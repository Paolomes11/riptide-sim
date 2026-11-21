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
// Sums of per-step deposits compared with each other or with a kinetic energy:
// each side is a sum of up to a few hundred doubles, so allow far more than
// their rounding but far less than any physical deposit.
constexpr double kSumAbsTolerance = 1e-9; // MeV
// G4NuclearStopping (11.4.1) caps its step loss at the pre-step kinetic energy,
// not at what ionisation left, and G4ParticleChangeForLoss then clamps the
// kinetic energy to 0 while still depositing the whole loss: on the last step
// of a stopping proton the track can deposit slightly more than it started
// with (1.9e-4 MeV seen at 14 MeV). The excess is bounded by one step of
// nuclear stopping; 10 keV is far above that and far below the MeV-scale error
// of tagging the wrong track.
constexpr double kNuclearStoppingSlack = 0.01; // MeV

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

const RunData& runKB0() {
    static const RunData run = loadRun(pathFromEnv("RIPTIDE_IT_RUN_KB0"));
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
    checkRecord(runKB0(), 1, pathFromEnv("RIPTIDE_IT_RUN_KB0"));
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

TEST_CASE("deposit and proton-track columns are consistent", "[integration]") {
    using riptide::schema::EventClass;
    const RunData& run = run2t();
    CHECK(run.record.at("schema_version").get<int>() == riptide::schema::kSchemaVersion);
    for (const auto& [event_id, row] : run.rows) {
        INFO("event_id " << event_id);
        const double Edep = doubleAt(row, "det_Edep");
        const double Evis = doubleAt(row, "det_Evis");
        CHECK(Edep >= 0.0);
        CHECK(Evis >= 0.0);
        // Birks quenching only ever reduces a step's deposit.
        CHECK(Evis <= Edep + kSumAbsTolerance);

        const auto event_class = intAt(row, "truth_event_class");
        if (static_cast<FirstInteraction>(intAt(row, "truth_first_interaction")) ==
            FirstInteraction::Unknown) {
            CHECK(event_class == static_cast<std::int32_t>(EventClass::NoInteraction));
            // Nothing interacted, so nothing was deposited: exactly zero.
            CHECK_THAT(Edep, Catch::Matchers::WithinAbs(0.0, 0.0));
        }
        if (intAt(row, "truth_valid_np") == 1) {
            const double p_Edep = doubleAt(row, "truth_p_Edep");
            const double p_Evis = doubleAt(row, "truth_p_Evis");
            CHECK(doubleAt(row, "truth_p_length") > 0.0);
            CHECK(p_Evis <= p_Edep + kSumAbsTolerance);
            CHECK(p_Edep <= doubleAt(row, "truth_Ep") + kNuclearStoppingSlack);
            CHECK(p_Edep <= Edep + kSumAbsTolerance);
            CHECK(event_class != static_cast<std::int32_t>(EventClass::NoInteraction));
        } else {
            CHECK(intAt(row, "truth_p_exits") == 0);
            CHECK(intAt(row, "truth_p_n_hadronic") == 0);
            CHECK(intAt(row, "truth_p_n_nuclear_recoil") == 0);
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

TEST_CASE("with kB = 0 the visible energy is the deposit", "[integration]") {
    const RunData& run = runKB0();
    REQUIRE_THAT(run.record.at("config").at("scintillator").at("birks_kB_mm_per_MeV").get<double>(),
                 Catch::Matchers::WithinAbs(0.0, 0.0));
    const auto& nominal = run1t().rows;
    REQUIRE(run.rows.size() == nominal.size());
    for (const auto& [event_id, row] : run.rows) {
        INFO("event_id " << event_id);
        // G4EmSaturation returns each step's deposit unchanged when kB = 0, and
        // both totals add the same steps in the same order.
        const double Edep = doubleAt(row, "det_Edep");
        CHECK_THAT(doubleAt(row, "det_Evis"), Catch::Matchers::WithinAbs(Edep, kSumAbsTolerance));
        CHECK_THAT(doubleAt(row, "truth_p_Evis"),
                   Catch::Matchers::WithinAbs(doubleAt(row, "truth_p_Edep"), kSumAbsTolerance));
        // kB only enters the visible energy, computed after the step: transport
        // and random numbers are those of the nominal 1-thread run.
        const auto it = nominal.find(event_id);
        REQUIRE(it != nominal.end());
        CHECK_THAT(Edep, Catch::Matchers::WithinULP(doubleAt(it->second, "det_Edep"), 0));
        CHECK(intAt(row, "truth_event_class") == intAt(it->second, "truth_event_class"));
    }
}
