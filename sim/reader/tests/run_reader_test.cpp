// Unit tests of the run reader on hand-written run.json files. Reading real
// ntuples is covered by the sim_integration test, which needs the data disk.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <nlohmann/json.hpp>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

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

TEST_CASE("run records of every readable schema version are accepted", "[run_reader]") {
    const auto dir = freshTempDir("riptide_run_reader_versions");
    for (std::int32_t v = riptide::schema::kOldestReadableSchemaVersion;
         v <= riptide::schema::kSchemaVersion; ++v) {
        INFO("schema_version " << v);
        json doc = validRecord();
        doc["schema_version"] = v;
        writeRecord(dir, doc);
        CHECK(riptide::sim::readRunRecord(dir) == doc);
    }
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
    SECTION("schema_version newer than this build") {
        doc["schema_version"] = riptide::schema::kSchemaVersion + 1;
    }
    SECTION("schema_version older than the oldest readable") {
        doc["schema_version"] = riptide::schema::kOldestReadableSchemaVersion - 1;
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
    const auto row = riptide::sim::eventRow(riptide::sim::EventMeta{kEventId, 1, 2}, truth,
                                            riptide::sim::DepositRecord{});

    CHECK(riptide::sim::intAt(row, "event_id") == kEventId);
    CHECK_THAT(riptide::sim::doubleAt(row, "truth_En"), Catch::Matchers::WithinULP(kEn, 0));
    CHECK_THROWS_AS(riptide::sim::intAt(row, "no_such_column"), std::out_of_range);
    CHECK_THROWS_AS(riptide::sim::intAt(row, "truth_En"), std::bad_variant_access);
    CHECK_THROWS_AS(riptide::sim::doubleAt(row, "event_id"), std::bad_variant_access);
}

namespace {

// Rows with the given event ids, handed out in that order.
class VectorRowSource final : public riptide::sim::RowSource {
public:
    VectorRowSource(std::string name, std::vector<std::int32_t> ids)
        : name_(std::move(name)), ids_(std::move(ids)) {}

    bool next(riptide::sim::EventRow& row) override {
        if (pos_ == ids_.size()) {
            return false;
        }
        row = riptide::sim::eventRow(riptide::sim::EventMeta{ids_[pos_], 0, 0},
                                     riptide::sim::TruthRecord{}, riptide::sim::DepositRecord{});
        ++pos_;
        return true;
    }

    [[nodiscard]] std::string name() const override {
        return name_;
    }

private:
    std::string name_;
    std::vector<std::int32_t> ids_;
    std::size_t pos_ = 0;
};

// The event ids in the order mergeByEventId visits them.
std::vector<std::int32_t> mergedIds(const std::vector<riptide::sim::RowSource*>& sources) {
    std::vector<std::int32_t> ids;
    const std::size_t n =
        riptide::sim::mergeByEventId(sources, [&ids](const riptide::sim::EventRow& row) {
            ids.push_back(riptide::sim::intAt(row, "event_id"));
        });
    CHECK(n == ids.size());
    return ids;
}

} // namespace

TEST_CASE("rows are merged in ascending event_id whatever the split", "[run_reader]") {
    const std::vector<std::int32_t> all{0, 1, 2, 3, 4, 5, 6, 7};
    SECTION("one source") {
        VectorRowSource a("a", all);
        CHECK(mergedIds({&a}) == all);
    }
    SECTION("interleaved chunks, as worker threads take them") {
        VectorRowSource a("a", {0, 1, 4, 5});
        VectorRowSource b("b", {2, 3, 6, 7});
        CHECK(mergedIds({&a, &b}) == all);
        VectorRowSource c("c", {2, 3, 6, 7});
        VectorRowSource d("d", {0, 1, 4, 5});
        CHECK(mergedIds({&c, &d}) == all);
    }
    SECTION("uneven split with an empty source") {
        VectorRowSource a("a", {});
        VectorRowSource b("b", {0, 2, 3, 4, 5, 6});
        VectorRowSource c("c", {1, 7});
        CHECK(mergedIds({&a, &b, &c}) == all);
    }
    SECTION("no sources") {
        CHECK(mergedIds({}).empty());
    }
}

TEST_CASE("an event id in two sources is visited twice, in source order", "[run_reader]") {
    VectorRowSource a("a", {0, 2});
    VectorRowSource b("b", {1, 2});
    CHECK(mergedIds({&a, &b}) == std::vector<std::int32_t>{0, 1, 2, 2});
}

TEST_CASE("a source not in ascending event_id is refused", "[run_reader]") {
    SECTION("descending") {
        VectorRowSource a("a", {0, 3, 2});
        CHECK_THROWS_AS(mergedIds({&a}), RunReadError);
    }
    SECTION("repeated within a source") {
        VectorRowSource a("a", {0, 1, 1});
        VectorRowSource b("b", {2});
        CHECK_THROWS_AS(mergedIds({&a, &b}), RunReadError);
    }
}
