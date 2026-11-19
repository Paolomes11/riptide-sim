// Reader side of the riptide_sim integration test (driven by
// run_integration.sh). It reads back two runs of the same events, written with
// 1 and 2 worker threads, through riptide::sim_reader (the Geant4 ROOT reader, so the simulation
// side keeps no ROOT dependency). These are consistency checks of the written files, not physics
// validation.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <variant>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/event_row.hpp"
#include "riptide/sim/run_reader.hpp"
#include "riptide/sim/run_record.hpp"

namespace {

using nlohmann::json;
using riptide::schema::ColumnType;
using riptide::schema::FirstInteraction;
using riptide::sim::columnSpec;
using riptide::sim::doubleAt;
using riptide::sim::intAt;
using riptide::sim::kEventColumns;
using riptide::sim::loadRun;
using riptide::sim::RunData;

// Direction cosines are normalised in double precision; this leaves room for a
// few rounding steps only.
constexpr double kUnitNormTolerance = 1e-12;
// The gun energy is written back unchanged, so the relative tolerance only
// absorbs the MeV conversion.
constexpr double kEnergyRelTolerance = 1e-12;

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

const RunData& run1t() {
    static const RunData run = loadRun(pathFromEnv("RIPTIDE_IT_RUN_1T"));
    return run;
}

const RunData& run2t() {
    static const RunData run = loadRun(pathFromEnv("RIPTIDE_IT_RUN_2T"));
    return run;
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
