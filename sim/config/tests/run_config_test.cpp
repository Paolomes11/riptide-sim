#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <string>

#include "riptide/sim/run_config.hpp"

namespace {

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;
using nlohmann::json;
using riptide::sim::ConfigError;
using riptide::sim::parseRunConfig;
using riptide::sim::SourceMode;

json validMono() {
    return json::parse(R"({
      "config_version": 2,
      "geometry": {"scint_size_mm": [100.0, 80.0, 60.0]},
      "source": {"mode": "mono", "En_MeV": 14.0,
                 "direction": [0.0, 0.0, 1.0], "position_mm": [0.0, 0.0, -100.0]},
      "physics_list": "QGSP_BIC_HP",
      "scintillator": {"birks_kB_mm_per_MeV": 0.155, "birks_source": "D18"},
      "proton_step_function": {"dRoverRange": 0.2, "finalRange_mm": 0.1}
    })");
}

// Layout of the phase 1b runs, still found in their run.json.
json validVersion1() {
    auto doc = validMono();
    doc["config_version"] = 1;
    doc.erase("scintillator");
    doc.erase("proton_step_function");
    return doc;
}

std::filesystem::path writeTemporary(const char* name, const std::string& text) {
    // A few bytes in the temporary directory: no run data involved.
    const auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream out(path);
    out << text;
    return path;
}

json validFlat() {
    auto doc = validMono();
    doc["source"] = json::parse(R"({"mode": "flat", "En_min_MeV": 10.0, "En_max_MeV": 20.0,
                                    "direction": [0.0, 0.0, 1.0],
                                    "position_mm": [0.0, 0.0, -100.0]})");
    return doc;
}

void requireError(const json& doc, const std::string& fragment) {
    REQUIRE_THROWS_MATCHES(parseRunConfig(doc), ConfigError,
                           Catch::Matchers::MessageMatches(ContainsSubstring(fragment)));
}

} // namespace

TEST_CASE("valid monoenergetic config is parsed in project units", "[config]") {
    const auto c = parseRunConfig(validMono());
    CHECK(c.config_version == riptide::sim::kConfigVersion);
    CHECK_THAT(c.geometry.scint_size.x, WithinULP(100.0, 0));
    CHECK_THAT(c.geometry.scint_size.y, WithinULP(80.0, 0));
    CHECK_THAT(c.geometry.scint_size.z, WithinULP(60.0, 0));
    CHECK(c.source.mode == SourceMode::Mono);
    CHECK_THAT(c.source.En_min, WithinULP(14.0, 0));
    CHECK_THAT(c.source.En_max, WithinULP(14.0, 0));
    CHECK_THAT(c.source.position.z, WithinULP(-100.0, 0));
    CHECK(c.physics_list == "QGSP_BIC_HP");
    REQUIRE(c.scintillator.has_value());
    CHECK_THAT(c.scintillator->birks_kB, WithinULP(0.155, 0));
    CHECK(c.scintillator->birks_source == "D18");
    REQUIRE(c.proton_step_function.has_value());
    CHECK_THAT(c.proton_step_function->dRoverRange, WithinULP(0.2, 0));
    CHECK_THAT(c.proton_step_function->finalRange, WithinULP(0.1, 0));
}

TEST_CASE("version 1 config is still readable, without scintillator physics", "[config]") {
    const auto c = parseRunConfig(validVersion1());
    CHECK(c.config_version == 1);
    CHECK_FALSE(c.scintillator.has_value());
    CHECK_FALSE(c.proton_step_function.has_value());
    CHECK(c.source.mode == SourceMode::Mono);
}

TEST_CASE("version 1 config refuses the version 2 keys", "[config]") {
    auto doc = validVersion1();
    doc["scintillator"] = validMono()["scintillator"];
    requireError(doc, "unknown key 'scintillator'");
}

TEST_CASE("zero Birks constant is allowed as the no-quenching oracle", "[config]") {
    auto doc = validMono();
    doc["scintillator"]["birks_kB_mm_per_MeV"] = 0.0;
    const auto c = parseRunConfig(doc);
    REQUIRE(c.scintillator.has_value());
    CHECK_THAT(c.scintillator->birks_kB, WithinAbs(0.0, 0.0));
}

TEST_CASE("valid flat config keeps both ends of the energy range", "[config]") {
    const auto c = parseRunConfig(validFlat());
    CHECK(c.source.mode == SourceMode::Flat);
    CHECK_THAT(c.source.En_min, WithinULP(10.0, 0));
    CHECK_THAT(c.source.En_max, WithinULP(20.0, 0));
}

TEST_CASE("direction within tolerance is normalised", "[config]") {
    auto doc = validMono();
    doc["source"]["direction"] = {0.0, 0.0, 1.0 + 5.0e-7};
    const auto d = parseRunConfig(doc).source.direction;
    CHECK_THAT(d.x * d.x + d.y * d.y + d.z * d.z, WithinAbs(1.0, 1.0e-15));
}

TEST_CASE("missing, unknown and mistyped keys are rejected", "[config]") {
    SECTION("missing top-level key") {
        auto doc = validMono();
        doc.erase("physics_list");
        requireError(doc, "missing key 'physics_list'");
    }
    SECTION("unknown top-level key") {
        auto doc = validMono();
        doc["seed"] = 1;
        requireError(doc, "unknown key 'seed'");
    }
    SECTION("misspelt source key") {
        auto doc = validMono();
        doc["source"]["En_Mev"] = 14.0;
        requireError(doc, "unknown key 'En_Mev'");
    }
    SECTION("flat key in mono mode") {
        auto doc = validMono();
        doc["source"]["En_max_MeV"] = 20.0;
        requireError(doc, "unknown key 'En_max_MeV'");
    }
    SECTION("unknown mode") {
        auto doc = validMono();
        doc["source"]["mode"] = "gauss";
        requireError(doc, "source.mode");
    }
    SECTION("number given as string") {
        auto doc = validMono();
        doc["source"]["En_MeV"] = "14";
        requireError(doc, "source.En_MeV");
    }
    SECTION("wrong vector length") {
        auto doc = validMono();
        doc["geometry"]["scint_size_mm"] = {100.0, 100.0};
        requireError(doc, "geometry.scint_size_mm");
    }
    SECTION("config version newer than the code") {
        auto doc = validMono();
        doc["config_version"] = riptide::sim::kConfigVersion + 1;
        requireError(doc, "config_version");
    }
    SECTION("config version older than any known") {
        auto doc = validMono();
        doc["config_version"] = riptide::sim::kOldestConfigVersion - 1;
        requireError(doc, "config_version");
    }
    SECTION("version 2 without scintillator") {
        auto doc = validMono();
        doc.erase("scintillator");
        requireError(doc, "missing key 'scintillator'");
    }
    SECTION("version 2 without step function") {
        auto doc = validMono();
        doc.erase("proton_step_function");
        requireError(doc, "missing key 'proton_step_function'");
    }
    SECTION("misspelt scintillator key") {
        auto doc = validMono();
        doc["scintillator"]["birks_kB"] = 0.155;
        requireError(doc, "unknown key 'birks_kB'");
    }
    SECTION("empty Birks source") {
        auto doc = validMono();
        doc["scintillator"]["birks_source"] = "";
        requireError(doc, "scintillator.birks_source");
    }
    SECTION("empty physics list") {
        auto doc = validMono();
        doc["physics_list"] = "";
        requireError(doc, "physics_list");
    }
}

TEST_CASE("physically invalid values are rejected", "[config]") {
    SECTION("non-positive edge") {
        auto doc = validMono();
        doc["geometry"]["scint_size_mm"] = {100.0, 0.0, 100.0};
        requireError(doc, "positive");
    }
    SECTION("zero energy") {
        auto doc = validMono();
        doc["source"]["En_MeV"] = 0.0;
        requireError(doc, "positive");
    }
    SECTION("energy above the evaluated-data limit") {
        auto doc = validMono();
        doc["source"]["En_MeV"] = 20.5;
        requireError(doc, "20 MeV");
    }
    SECTION("20 MeV itself is allowed") {
        auto doc = validMono();
        doc["source"]["En_MeV"] = 20.0;
        CHECK_NOTHROW(parseRunConfig(doc));
    }
    SECTION("empty flat range") {
        auto doc = validFlat();
        doc["source"]["En_min_MeV"] = 20.0;
        requireError(doc, "smaller");
    }
    SECTION("non-unit direction") {
        auto doc = validMono();
        doc["source"]["direction"] = {0.0, 0.0, 2.0};
        requireError(doc, "unit vector");
    }
    SECTION("negative Birks constant") {
        auto doc = validMono();
        doc["scintillator"]["birks_kB_mm_per_MeV"] = -0.1;
        requireError(doc, "negative");
    }
    SECTION("step function ratio outside (0, 1]") {
        for (const double ratio : {0.0, -0.1, 1.5}) {
            auto doc = validMono();
            doc["proton_step_function"]["dRoverRange"] = ratio;
            requireError(doc, "proton_step_function.dRoverRange");
        }
    }
    SECTION("step function ratio 1 is allowed") {
        auto doc = validMono();
        doc["proton_step_function"]["dRoverRange"] = 1.0;
        CHECK_NOTHROW(parseRunConfig(doc));
    }
    SECTION("non-positive final range") {
        auto doc = validMono();
        doc["proton_step_function"]["finalRange_mm"] = 0.0;
        requireError(doc, "proton_step_function.finalRange_mm");
    }
    SECTION("source inside the scintillator") {
        auto doc = validMono();
        doc["source"]["position_mm"] = {0.0, 0.0, -20.0};
        requireError(doc, "outside");
    }
}

TEST_CASE("toJson round-trips through parseRunConfig", "[config]") {
    for (const auto& doc : {validMono(), validFlat(), validVersion1()}) {
        const auto c = parseRunConfig(doc);
        const auto back = parseRunConfig(riptide::sim::toJson(c));
        // Same document both ways, version 1 included: nothing is lost or added.
        CHECK(riptide::sim::toJson(back) == riptide::sim::toJson(c));
        CHECK(back.config_version == c.config_version);
        CHECK(back.scintillator.has_value() == c.scintillator.has_value());
        CHECK(back.proton_step_function.has_value() == c.proton_step_function.has_value());
        CHECK(back.source.mode == c.source.mode);
        CHECK_THAT(back.source.En_min, WithinULP(c.source.En_min, 0));
        CHECK_THAT(back.source.En_max, WithinULP(c.source.En_max, 0));
        CHECK_THAT(back.geometry.scint_size.y, WithinULP(c.geometry.scint_size.y, 0));
        CHECK_THAT(back.source.direction.z, WithinULP(c.source.direction.z, 0));
        CHECK_THAT(back.source.position.z, WithinULP(c.source.position.z, 0));
        CHECK(back.physics_list == c.physics_list);
    }
}

TEST_CASE("shipped config files are valid", "[config]") {
    const std::filesystem::path dir{RIPTIDE_CONFIG_DIR};
    const auto mono = riptide::sim::loadRunConfig(dir / "mono_14MeV.json");
    CHECK(mono.source.mode == SourceMode::Mono);
    const auto flat = riptide::sim::loadRunConfig(dir / "flat_10_20MeV.json");
    CHECK(flat.source.mode == SourceMode::Flat);
    // The edge run covers exactly the HP/CHIPS elastic overlap up to 20 MeV.
    const auto edge = riptide::sim::loadRunConfig(dir / "flat_19p5_20MeV.json");
    CHECK(edge.source.mode == SourceMode::Flat);
    CHECK_THAT(edge.source.En_min, WithinULP(19.5, 0));
    CHECK_THAT(edge.source.En_max, WithinULP(20.0, 0));
    // Every shipped run uses the Birks constant of D18 and Geant4's default
    // step function for muons and hadrons (G4EmExtraParameters, 11.4.1).
    for (const auto& c : {mono, flat, edge}) {
        CHECK(c.config_version == riptide::sim::kConfigVersion);
        REQUIRE(c.scintillator.has_value());
        CHECK_THAT(c.scintillator->birks_kB, WithinULP(0.155, 0));
        CHECK(c.scintillator->birks_source == "D18");
        REQUIRE(c.proton_step_function.has_value());
        CHECK_THAT(c.proton_step_function->dRoverRange, WithinULP(0.2, 0));
        CHECK_THAT(c.proton_step_function->finalRange, WithinULP(0.1, 0));
    }
}

TEST_CASE("a new run refuses an older config version", "[config]") {
    const auto path = writeTemporary("riptide_config_test_v1.json", validVersion1().dump());
    CHECK_THROWS_MATCHES(riptide::sim::loadRunConfig(path), ConfigError,
                         Catch::Matchers::MessageMatches(ContainsSubstring("only readable")));
    std::filesystem::remove(path);
}

TEST_CASE("unreadable or malformed files raise ConfigError", "[config]") {
    CHECK_THROWS_AS(riptide::sim::loadRunConfig("/nonexistent/riptide.json"), ConfigError);

    const auto path = writeTemporary("riptide_config_test_bad.json", "{\"config_version\": 2,");
    CHECK_THROWS_MATCHES(riptide::sim::loadRunConfig(path), ConfigError,
                         Catch::Matchers::MessageMatches(ContainsSubstring("invalid JSON")));
    std::filesystem::remove(path);
}

TEST_CASE("world contains the scintillator and the gun", "[config]") {
    const auto c = parseRunConfig(validMono()); // box 100 x 80 x 60, gun at z = -100
    const double half = riptide::sim::worldHalfSize(c);
    CHECK(half > 100.0);
    CHECK_THAT(half, WithinULP(riptide::sim::kWorldScale * 100.0 + riptide::sim::kWorldMargin, 0));
}
