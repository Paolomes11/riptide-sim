// Reader side of the riptide_sim integration test (driven by
// run_integration.sh). It reads back two runs of the same events, written with
// 1 and 2 worker threads, through the Geant4 ROOT reader: the simulation side
// keeps no ROOT dependency. These are consistency checks of the written files,
// not physics validation.

#include <G4RootAnalysisReader.hh>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/event_row.hpp"
#include "riptide/sim/run_record.hpp"

namespace {

using nlohmann::json;
using riptide::schema::ColumnType;
using riptide::schema::FirstInteraction;
using riptide::sim::columnIndex;
using riptide::sim::columnSpec;
using riptide::sim::ColumnValue;
using riptide::sim::EventRow;
using riptide::sim::kEventColumns;

// Direction cosines are normalised in double precision; this leaves room for a
// few rounding steps only.
constexpr double kUnitNormTolerance = 1e-12;
// The gun energy is written back unchanged, so the relative tolerance only
// absorbs the MeV conversion.
constexpr double kEnergyRelTolerance = 1e-12;
// G4AnalysisReader's "no ntuple" id (kInvalidId in G4AnalysisUtilities.hh).
constexpr G4int kInvalidNtupleId = -1;

struct RunData {
    json record;
    std::map<std::int32_t, EventRow> rows;
    std::size_t row_count = 0; // rows read, duplicates included
};

std::filesystem::path pathFromEnv(const char* name) {
    const char* value = std::getenv(name); // NOLINT(concurrency-mt-unsafe): single-threaded
    if (value == nullptr || *value == '\0') {
        throw std::runtime_error(std::string(name) + " is not set; run through run_integration.sh");
    }
    return value;
}

std::size_t expectedEvents() {
    return static_cast<std::size_t>(std::stoul(pathFromEnv("RIPTIDE_IT_EVENTS").string()));
}

// Reads one per-thread file into run.rows. Every schema column is bound by
// name, so a column missing from the file fails here.
void readEventFile(const std::filesystem::path& file, RunData& run) {
    auto* reader = G4RootAnalysisReader::Instance();
    const G4int id = reader->GetNtuple(G4String(riptide::schema::kEventTreeName), file.string());
    if (id == kInvalidNtupleId) {
        throw std::runtime_error("cannot read ntuple from " + file.string());
    }
    std::array<G4int, kEventColumns> ints{};
    std::array<G4double, kEventColumns> doubles{};
    for (std::size_t i = 0; i < kEventColumns; ++i) {
        const auto& spec = columnSpec(i);
        const G4String name(spec.name);
        const bool bound = spec.type == ColumnType::Int32
                               ? reader->SetNtupleIColumn(id, name, ints.at(i))
                               : reader->SetNtupleDColumn(id, name, doubles.at(i));
        if (!bound) {
            throw std::runtime_error("cannot bind column " + name + " in " + file.string());
        }
    }
    const std::size_t event_id_col = columnIndex("event_id");
    while (reader->GetNtupleRow(id)) {
        EventRow row;
        for (std::size_t i = 0; i < kEventColumns; ++i) {
            row.at(i) = columnSpec(i).type == ColumnType::Int32 ? ColumnValue{ints.at(i)}
                                                                : ColumnValue{doubles.at(i)};
        }
        ++run.row_count;
        run.rows.emplace(ints.at(event_id_col), row);
    }
}

RunData loadRun(const std::filesystem::path& dir) {
    std::ifstream in(dir / riptide::sim::kRunRecordFileName);
    if (!in) {
        throw std::runtime_error("no run record in " + dir.string());
    }
    RunData run;
    run.record = json::parse(in);
    if (run.record.at("schema_version").get<std::int32_t>() != riptide::schema::kSchemaVersion) {
        throw std::runtime_error("unknown schema_version in " + dir.string());
    }
    for (const auto& name : run.record.at("event_files")) {
        readEventFile(dir / name.get<std::string>(), run);
    }
    G4RootAnalysisReader::Instance()->CloseFiles();
    return run;
}

const RunData& run1t() {
    static const RunData run = loadRun(pathFromEnv("RIPTIDE_IT_RUN_1T"));
    return run;
}

const RunData& run2t() {
    static const RunData run = loadRun(pathFromEnv("RIPTIDE_IT_RUN_2T"));
    return run;
}

std::int32_t intAt(const EventRow& row, std::string_view name) {
    return std::get<std::int32_t>(row.at(columnIndex(name)));
}

double doubleAt(const EventRow& row, std::string_view name) {
    return std::get<double>(row.at(columnIndex(name)));
}

void checkRecord(const RunData& run, int threads, const std::filesystem::path& dir) {
    CHECK(run.record.at("events").get<std::size_t>() == expectedEvents());
    CHECK(run.record.at("threads").get<int>() == threads);
    CHECK(run.record.at("event_files") == json(riptide::sim::eventFileNames(threads)));
    for (const auto& name : run.record.at("event_files")) {
        CHECK(std::filesystem::is_regular_file(dir / name.get<std::string>()));
    }
}

} // namespace

TEST_CASE("run records match the requested runs", "[integration]") {
    checkRecord(run1t(), 1, pathFromEnv("RIPTIDE_IT_RUN_1T"));
    checkRecord(run2t(), 2, pathFromEnv("RIPTIDE_IT_RUN_2T"));
}

TEST_CASE("event ids are unique and cover 0..N-1", "[integration]") {
    const std::size_t n = expectedEvents();
    for (const RunData* run : {&run1t(), &run2t()}) {
        CHECK(run->row_count == n);
        REQUIRE(run->rows.size() == n);
        // Keys are sorted and unique, so first == 0 and last == N-1 means 0..N-1.
        CHECK(run->rows.begin()->first == 0);
        CHECK(static_cast<std::size_t>(run->rows.rbegin()->first) == n - 1);
    }
}

TEST_CASE("truth columns are consistent", "[integration]") {
    const RunData& run = run2t();
    const double En_mono = run.record.at("config").at("source").at("En_MeV").get<double>();
    for (const auto& [event_id, row] : run.rows) {
        INFO("event_id " << event_id);
        const double En = doubleAt(row, "truth_En");
        CHECK_THAT(En, Catch::Matchers::WithinRel(En_mono, kEnergyRelTolerance));

        const double norm =
            std::hypot(doubleAt(row, "truth_omega_n_x"), doubleAt(row, "truth_omega_n_y"),
                       doubleAt(row, "truth_omega_n_z"));
        CHECK_THAT(norm, Catch::Matchers::WithinAbs(1.0, kUnitNormTolerance));

        const auto first = static_cast<FirstInteraction>(intAt(row, "truth_first_interaction"));
        const bool valid_np = intAt(row, "truth_valid_np") == 1;
        if (valid_np) {
            const double Ep = doubleAt(row, "truth_Ep");
            CHECK(Ep > 0.0);
            CHECK(Ep <= En);
            // Only a first interaction on H keeps the neutron on the beam axis;
            // after an earlier n-C scatter theta_p can exceed pi/2.
            if (first == FirstInteraction::NpElastic) {
                CHECK(doubleAt(row, "truth_theta_p") <= riptide::units::pi / 2.0);
            }
        }
        if (first == FirstInteraction::Unknown) {
            CHECK_FALSE(valid_np);
            CHECK(intAt(row, "truth_n_np_scatters") == 0);
        }
    }
}

TEST_CASE("the same events come out with 1 and 2 threads", "[integration]") {
    const auto& rows1 = run1t().rows;
    const auto& rows2 = run2t().rows;
    REQUIRE(rows1.size() == rows2.size());
    for (const auto& [event_id, row1] : rows1) {
        INFO("event_id " << event_id);
        const auto it = rows2.find(event_id);
        REQUIRE(it != rows2.end());
        for (std::size_t i = 0; i < kEventColumns; ++i) {
            INFO("column " << columnSpec(i).name);
            if (columnSpec(i).type == ColumnType::Int32) {
                CHECK(std::get<std::int32_t>(row1.at(i)) ==
                      std::get<std::int32_t>(it->second.at(i)));
            } else {
                CHECK_THAT(std::get<double>(row1.at(i)),
                           Catch::Matchers::WithinULP(std::get<double>(it->second.at(i)), 0));
            }
        }
    }
}
