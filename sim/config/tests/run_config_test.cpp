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
      "config_version": 1,
      "geometry": {"scint_size_mm": [100.0, 80.0, 60.0]},
      "source": {"mode": "mono", "En_MeV": 14.0,
                 "direction": [0.0, 0.0, 1.0], "position_mm": [0.0, 0.0, -100.0]},
      "physics_list": "QGSP_BIC_HP"
    })");
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
    SECTION("wrong config version") {
        auto doc = validMono();
        doc["config_version"] = 2;
        requireError(doc, "config_version");
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
    SECTION("source inside the scintillator") {
        auto doc = validMono();
        doc["source"]["position_mm"] = {0.0, 0.0, -20.0};
        requireError(doc, "outside");
    }
}

TEST_CASE("toJson round-trips through parseRunConfig", "[config]") {
    for (const auto& doc : {validMono(), validFlat()}) {
        const auto c = parseRunConfig(doc);
        const auto back = parseRunConfig(riptide::sim::toJson(c));
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
}

TEST_CASE("unreadable or malformed files raise ConfigError", "[config]") {
    CHECK_THROWS_AS(riptide::sim::loadRunConfig("/nonexistent/riptide.json"), ConfigError);

    // A few bytes in the temporary directory: no run data involved.
    const auto path = std::filesystem::temp_directory_path() / "riptide_config_test_bad.json";
    {
        std::ofstream out(path);
        out << "{\"config_version\": 1,";
    }
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
