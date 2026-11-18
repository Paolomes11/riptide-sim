// Unit tests of the run reader on hand-written run.json files. Reading real
// ntuples is covered by the sim_integration test, which needs the data disk.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <nlohmann/json.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <variant>

#include "riptide/common/event_schema.hpp"
#include "riptide/sim/event_row.hpp"
#include "riptide/sim/run_reader.hpp"
#include "riptide/sim/run_record.hpp"
#include "riptide/sim/truth_record.hpp"

namespace {

using nlohmann::json;
using riptide::sim::RunReadError;

// A fresh, empty directory under the system temporary directory (a few bytes
// are written there; no run data involved).
std::filesystem::path freshTempDir(const std::string& name) {
    const auto dir = std::filesystem::temp_directory_path() / name;
    std::filesystem::remove_all(dir);
    std::filesystem::create_directory(dir);
    return dir;
}

// The keys readRunRecord checks, with the values this build writes.
json validRecord() {
    json doc;
    doc["run_record_version"] = riptide::sim::kRunRecordVersion;
    doc["schema_version"] = riptide::schema::kSchemaVersion;
    doc["event_files"] = json::array({"events_t0.root"});
    return doc;
}

void writeRecord(const std::filesystem::path& dir, const json& doc) {
    std::ofstream(dir / riptide::sim::kRunRecordFileName) << doc.dump();
}

} // namespace

TEST_CASE("a valid run record is read back", "[run_reader]") {
    const auto dir = freshTempDir("riptide_run_reader_valid");
    writeRecord(dir, validRecord());
    CHECK(riptide::sim::readRunRecord(dir) == validRecord());
}

TEST_CASE("a run without run.json is refused", "[run_reader]") {
    const auto dir = freshTempDir("riptide_run_reader_missing");
    CHECK_THROWS_AS(riptide::sim::readRunRecord(dir), RunReadError);
}

TEST_CASE("run records of another layout are refused", "[run_reader]") {
    const auto dir = freshTempDir("riptide_run_reader_layout");
    json doc = validRecord();

    SECTION("run_record_version") {
        doc["run_record_version"] = riptide::sim::kRunRecordVersion + 1;
    }
    SECTION("schema_version") {
        doc["schema_version"] = riptide::schema::kSchemaVersion + 1;
    }
    SECTION("event_files missing") {
        doc.erase("event_files");
    }
    SECTION("event_files not a list") {
        doc["event_files"] = "events_t0.root";
    }
    SECTION("event_files with a non-string") {
        doc["event_files"] = json::array({0});
    }

    writeRecord(dir, doc);
    CHECK_THROWS_AS(riptide::sim::readRunRecord(dir), RunReadError);
}

TEST_CASE("malformed run.json is refused", "[run_reader]") {
    const auto dir = freshTempDir("riptide_run_reader_malformed");
    std::ofstream(dir / riptide::sim::kRunRecordFileName) << "{ not json";
    CHECK_THROWS_AS(riptide::sim::readRunRecord(dir), RunReadError);
}

TEST_CASE("a listed event file that does not exist is an error", "[run_reader]") {
    const auto dir = freshTempDir("riptide_run_reader_no_file");
    writeRecord(dir, validRecord());
    const json record = riptide::sim::readRunRecord(dir);
    std::size_t visited = 0;
    CHECK_THROWS_AS(
        riptide::sim::forEachEventRow(dir, record, [&visited](const auto&) { ++visited; }),
        RunReadError);
    CHECK(visited == 0);
}

TEST_CASE("columns are read by name and type", "[run_reader]") {
    constexpr std::int32_t kEventId = 42;
    constexpr double kEn = 14.0; // MeV
    riptide::sim::TruthRecord truth;
    truth.En = kEn;
    const auto row = riptide::sim::eventRow(riptide::sim::EventMeta{kEventId, 1, 2}, truth);

    CHECK(riptide::sim::intAt(row, "event_id") == kEventId);
    CHECK_THAT(riptide::sim::doubleAt(row, "truth_En"), Catch::Matchers::WithinULP(kEn, 0));
    CHECK_THROWS_AS(riptide::sim::intAt(row, "no_such_column"), std::out_of_range);
    CHECK_THROWS_AS(riptide::sim::intAt(row, "truth_En"), std::bad_variant_access);
    CHECK_THROWS_AS(riptide::sim::doubleAt(row, "event_id"), std::bad_variant_access);
}
