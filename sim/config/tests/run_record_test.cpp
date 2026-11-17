#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/version.hpp"
#include "riptide/sim/run_config.hpp"
#include "riptide/sim/run_record.hpp"

namespace {

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinULP;
using nlohmann::json;
using riptide::sim::RunProvenance;
using riptide::sim::RunRecordError;

riptide::sim::RunConfig shippedConfig(const char* name) {
    return riptide::sim::loadRunConfig(std::filesystem::path{RIPTIDE_CONFIG_DIR} / name);
}

RunProvenance sampleRun() {
    RunProvenance run;
    run.run_id = "unit_test";
    run.seed = 20260907;
    run.events = 1000;
    run.threads = 2;
    run.geant4_version = "11.4.1";
    run.started_utc = "2026-09-07T12:34:56Z";
    run.elapsed_s = 1.5;
    return run;
}

// A fresh, empty directory under the system temporary directory (a few bytes
// are written there; no run data involved).
std::filesystem::path freshTempDir(const std::string& name) {
    const auto dir = std::filesystem::temp_directory_path() / name;
    std::filesystem::remove_all(dir);
    std::filesystem::create_directory(dir);
    return dir;
}

} // namespace

TEST_CASE("formatUtc writes ISO 8601 in UTC", "[run_record]") {
    CHECK(riptide::sim::formatUtc(0) == "1970-01-01T00:00:00Z");
    constexpr std::time_t kKnownTime = 1788784496; // 2026-09-07 12:34:56 UTC
    CHECK(riptide::sim::formatUtc(kKnownTime) == "2026-09-07T12:34:56Z");
}

TEST_CASE("event file names follow the Geant4 per-thread rule", "[run_record]") {
    CHECK(riptide::sim::eventFileNames(1) == std::vector<std::string>{"events_t0.root"});
    CHECK(riptide::sim::eventFileNames(3) ==
          std::vector<std::string>{"events_t0.root", "events_t1.root", "events_t2.root"});
}

TEST_CASE("run record holds provenance with the expected types", "[run_record]") {
    const auto config = shippedConfig("mono_14MeV.json");
    const auto doc = riptide::sim::runRecordJson(config, sampleRun());

    const std::vector<std::string> keys{"run_record_version",
                                        "run_id",
                                        "project_version",
                                        "git_describe",
                                        "geant4_version",
                                        "schema_version",
                                        "config_version",
                                        "seed",
                                        "events",
                                        "threads",
                                        "started_utc",
                                        "elapsed_s",
                                        "event_files",
                                        "config"};
    CHECK(doc.size() == keys.size());
    for (const auto& key : keys) {
        INFO("key " << key);
        CHECK(doc.contains(key));
    }

    CHECK(doc["run_record_version"].get<int>() == riptide::sim::kRunRecordVersion);
    CHECK(doc["run_id"].get<std::string>() == "unit_test");
    CHECK(doc["project_version"].get<std::string>() == riptide::kProjectVersion);
    CHECK(doc["git_describe"].get<std::string>() == riptide::kGitDescribe);
    CHECK_FALSE(doc["git_describe"].get<std::string>().empty());
    CHECK(doc["geant4_version"].get<std::string>() == "11.4.1");
    CHECK(doc["schema_version"].get<int>() == riptide::schema::kSchemaVersion);
    CHECK(doc["config_version"].get<int>() == riptide::sim::kConfigVersion);
    CHECK(doc["seed"].is_number_unsigned());
    CHECK(doc["seed"].get<std::uint64_t>() == 20260907U);
    CHECK(doc["events"].get<std::uint64_t>() == 1000U);
    CHECK(doc["threads"].get<int>() == 2);
    CHECK(doc["started_utc"].get<std::string>() == "2026-09-07T12:34:56Z");
    CHECK_THAT(doc["elapsed_s"].get<double>(), WithinULP(1.5, 0));
    CHECK(doc["event_files"] == json::array({"events_t0.root", "events_t1.root"}));
}

TEST_CASE("run record config parses back to the run config", "[run_record]") {
    for (const char* name : {"mono_14MeV.json", "flat_10_20MeV.json"}) {
        INFO("config " << name);
        const auto config = shippedConfig(name);
        const auto doc = riptide::sim::runRecordJson(config, sampleRun());
        const auto back = riptide::sim::parseRunConfig(doc["config"]);
        CHECK(back.source.mode == config.source.mode);
        CHECK_THAT(back.source.En_min, WithinULP(config.source.En_min, 0));
        CHECK_THAT(back.source.En_max, WithinULP(config.source.En_max, 0));
        CHECK_THAT(back.geometry.scint_size.x, WithinULP(config.geometry.scint_size.x, 0));
        CHECK(back.physics_list == config.physics_list);
        // Same document both ways: nothing is lost or added by the round trip.
        CHECK(riptide::sim::toJson(back) == doc["config"]);
    }
}

TEST_CASE("run record is written atomically and never overwritten", "[run_record]") {
    const auto dir = freshTempDir("riptide_run_record_test");
    const auto doc = riptide::sim::runRecordJson(shippedConfig("mono_14MeV.json"), sampleRun());

    riptide::sim::writeRunRecord(dir, doc);
    const auto path = dir / riptide::sim::kRunRecordFileName;
    REQUIRE(std::filesystem::exists(path));
    CHECK_FALSE(std::filesystem::exists(dir / "run.json.tmp"));
    {
        std::ifstream in(path);
        CHECK(json::parse(in) == doc);
    }

    CHECK_THROWS_MATCHES(riptide::sim::writeRunRecord(dir, doc), RunRecordError,
                         Catch::Matchers::MessageMatches(ContainsSubstring("already exists")));
    std::filesystem::remove_all(dir);
}

TEST_CASE("run record write fails in a missing directory", "[run_record]") {
    const auto dir = std::filesystem::temp_directory_path() / "riptide_run_record_missing";
    std::filesystem::remove_all(dir);
    CHECK_THROWS_AS(riptide::sim::writeRunRecord(dir, json::object()), RunRecordError);
    CHECK_FALSE(std::filesystem::exists(dir));
}
