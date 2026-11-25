// Calibrates the Case B engine with pseudo-experiments drawn from the toy
// model (D43), so the port can be compared with the toy before any Geant4
// data: the checklist on the (N, M) pairs of the toy plus a repeat at
// N = 1000, the contraction of the error with N, and the robustness to a vMF
// prior on the direction and to a wider E_n prior. Experiments run on
// several threads (D44); each has its own engines and output slot, so the
// numbers do not depend on --threads. It writes case_b_calibration.json under
// the stats directory of the data root; the file holds no wall time, so two
// runs with the same settings give the same bytes. Timings go to stdout only.

#include <lyra/lyra.hpp>
#include <nlohmann/json.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "riptide/apps/case_a_geant4.hpp"
#include "riptide/common/data_paths.hpp"
#include "riptide/common/direction.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/units.hpp"
#include "riptide/common/version.hpp"
#include "riptide/io/analysis_record.hpp"
#include "riptide/stats/calibration_b.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_b_posterior.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"

namespace {

namespace fs = std::filesystem;
namespace st = riptide::stats;
using Clock = std::chrono::steady_clock;
using nlohmann::json;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
constexpr std::uint64_t kProjectSeed = 20260907;
constexpr st::CaseBResolution kToyResolution{0.10, 0.08}; // MeV, rad
constexpr unsigned kDefaultThreads = 4;
constexpr int kCalibrationVersion = 1;
constexpr const char* kCalibrationName = "case_b_calibration.json";
// Generous upper bound on the output: one small JSON file.
constexpr std::uintmax_t kOutputBytes = std::uintmax_t{1} << 20U; // 1 MiB

// The protocol of the toy checklist, which the numbers are compared with.
struct SizedRow {
    std::size_t events;      // N per experiment
    std::size_t experiments; // M
};
const std::vector<SizedRow> kChecklistRows{
    {10, 400}, {30, 300}, {100, 200}, {300, 100}, {1000, 40}};
// The repeat at N = 1000 on disjoint ids, for the coverage the M = 40 row
// cannot resolve.
constexpr SizedRow kRepeatRow{1000, 240};
// The contraction slope is fitted from this N up, where the posterior is
// close to Gaussian.
constexpr std::size_t kSlopeMinEvents = 30;
// The robustness rows: the checklist rows from this N up, on their very
// experiments.
constexpr std::size_t kPriorMinEvents = 100;
constexpr double kVmfKappa = 1.0;
constexpr double kVmfAxisAngle = 0.5 * riptide::units::pi; // rad: steepest log-prior
constexpr double kWideEnMax = 30.0;                        // MeV: E_n prior on [10, 30]

// Disjoint event-id sections, one per row, so no two draws share an engine.
constexpr std::uint64_t kSectionChecklist = 16;
constexpr std::uint64_t kSectionRepeat = 24;
constexpr std::uint64_t kSectionPriorAxis = 32;

struct Options {
    std::string data_dir;
    std::string calibration_id;
    std::uint64_t seed = kProjectSeed;
    st::CaseBResolution resolution = kToyResolution;
    unsigned threads = kDefaultThreads;
    bool pilot = false;
};

class Stopwatch {
public:
    explicit Stopwatch(std::string name) : name_(std::move(name)) {}
    void report(std::size_t experiments) const {
        const std::chrono::duration<double> elapsed = Clock::now() - start_;
        std::cout << "  " << name_ << ": " << elapsed.count() << " s ("
                  << elapsed.count() / static_cast<double>(experiments) << " s per experiment)\n"
                  << std::flush;
    }

private:
    std::string name_;
    Clock::time_point start_ = Clock::now();
};

json domainJson(const st::EnergyDomain& d) {
    return {{"En_min_MeV", d.En_min}, {"En_max_MeV", d.En_max}, {"points", d.points}};
}

json ensembleJson(const st::CaseBEnsemble& e) {
    return {{"entries", e.entries},
            {"error_mean_rad", e.error_mean},
            {"error_mean_error_rad", e.error_mean_error},
            {"error_rms_rad", e.error_rms},
            {"sigma_mean_rad", e.sigma_mean},
            {"sigma_mean_error_rad", e.sigma_mean_error},
            {"pull_rms", e.pull_rms},
            {"scaled_pull_median", e.scaled_pull_median},
            {"scaled_pull_tail_fraction", e.scaled_pull_tail_fraction},
            {"truth_in_cap_fraction", e.truth_in_cap_fraction},
            {"coverage", e.coverage}};
}

json shiftJson(const st::DirectionPriorShift& s) {
    return {{"shift_mean", s.shift_mean},
            {"shift_error", s.shift_error},
            {"shift_max", s.shift_max},
            {"vector_mean", s.vector_mean},
            {"vector_error", s.vector_error}};
}

// The models and grids of the calibration: the reference model (E_n prior
// flat on the truth domain) and the alternative with the wider E_n prior.
struct Setup {
    st::EnergyDomain domain;
    st::EnergyDomain wide = st::sameStepDomain(domain, kWideEnMax);
    st::CaseBModel model;
    st::CaseBModel wide_model;
    std::vector<riptide::geometry::Direction> coarse_grid = st::sphereGrid(st::kSphereGridPixels);

    static st::CaseBModel modelOn(const st::EnergyDomain& d, const st::CaseBResolution& r) {
        const st::Grid En = st::energyGrid(d);
        return {En, st::recoilAngleGrid(st::kTrackAngleGridPoints),
                st::Grid{0.0, riptide::units::pi, st::kTrackObservedGridPoints},
                st::flatEnergyLogPrior(En), r};
    }
    explicit Setup(const st::CaseBResolution& r)
        : model(modelOn(domain, r)), wide_model(modelOn(wide, r)) {}
};

struct RowResult {
    st::CaseBRowSpec spec;
    std::vector<st::CaseBExperimentResult> results;
    st::CaseBEnsemble ensemble;
};

RowResult runRow(const Setup& s, const Options& o, const SizedRow& row, std::uint64_t section,
                 const std::string& name) {
    const st::CaseBRowSpec spec{row.events, riptide::apps::pilotScaled(row.experiments, o.pilot),
                                section};
    const Stopwatch watch(name + " N = " + std::to_string(spec.events) +
                          ", M = " + std::to_string(spec.experiments));
    RowResult out{
        spec,
        st::runCaseBRow(s.model, s.domain, spec, o.resolution, o.seed, s.coarse_grid, o.threads),
        {}};
    out.ensemble = st::summariseCaseBEnsemble(out.results);
    watch.report(spec.experiments);
    return out;
}

json rowJson(const RowResult& r) {
    const double root_n = std::sqrt(static_cast<double>(r.spec.events));
    return {{"events", r.spec.events},
            {"experiments", r.spec.experiments},
            {"section", r.spec.section},
            {"ensemble", ensembleJson(r.ensemble)},
            {"error_rms_sqrt_N_rad", r.ensemble.error_rms * root_n}};
}

// log-log slope of the rms error over the rows from kSlopeMinEvents up.
double slopeOf(const std::vector<const RowResult*>& rows) {
    std::vector<double> events;
    std::vector<double> rms;
    for (const RowResult* r : rows) {
        if (r->spec.events >= kSlopeMinEvents) {
            events.push_back(static_cast<double>(r->spec.events));
            rms.push_back(r->ensemble.error_rms);
        }
    }
    return st::logLogSlope(events, rms);
}

json priorRows(const Setup& s, const Options& o, const std::vector<RowResult>& checklist) {
    json rows = json::array();
    for (std::size_t r = 0; r < checklist.size(); ++r) {
        const RowResult& ref = checklist[r];
        if (ref.spec.events < kPriorMinEvents) {
            continue;
        }
        const Stopwatch watch("priors N = " + std::to_string(ref.spec.events) +
                              ", M = " + std::to_string(ref.spec.experiments));
        const st::CaseBPriorSpec prior{kVmfKappa, kVmfAxisAngle, kSectionPriorAxis + r};
        const st::CaseBPriorRow alt =
            st::runCaseBPriorRow(s.model, s.wide_model, s.domain, ref.spec, prior, o.resolution,
                                 o.seed, s.coarse_grid, o.threads);
        rows.push_back({{"events", ref.spec.events},
                        {"experiments", ref.spec.experiments},
                        {"vmf",
                         {{"kappa", prior.kappa},
                          {"axis_angle_rad", prior.axis_angle},
                          {"axis_section", prior.axis_section},
                          {"shift", shiftJson(st::directionPriorShift(ref.results, alt.vmf))},
                          {"ensemble", ensembleJson(st::summariseCaseBEnsemble(alt.vmf))}}},
                        {"energy",
                         {{"prior_domain", domainJson(s.wide)},
                          {"shift", shiftJson(st::directionPriorShift(ref.results, alt.energy))},
                          {"ensemble", ensembleJson(st::summariseCaseBEnsemble(alt.energy))}}}});
        watch.report(2 * ref.spec.experiments);
    }
    return rows;
}

int calibrate(const Options& o) {
    namespace data = riptide::data;
    const std::optional<std::string> cli_root =
        o.data_dir.empty() ? std::nullopt : std::optional<std::string>{o.data_dir};
    const fs::path root = data::resolveDataRoot(cli_root);
    const fs::path out_dir = data::statsDir(root, o.calibration_id);
    if (fs::exists(out_dir)) {
        throw data::DataPathError("calibration directory already exists: " + out_dir.string());
    }
    data::requireWritableRoot(root, kOutputBytes, data::kFreeSpaceMargin, false);

    std::cout << "riptide_case_b_calibrate " << riptide::kProjectVersion << " ("
              << riptide::kGitDescribe << "), " << o.threads << " threads"
              << (o.pilot ? ", pilot" : "") << '\n';
    const Setup setup(o.resolution);

    std::vector<RowResult> checklist;
    for (std::size_t r = 0; r < kChecklistRows.size(); ++r) {
        checklist.push_back(
            runRow(setup, o, kChecklistRows[r], kSectionChecklist + r, "checklist"));
    }
    const RowResult repeat = runRow(setup, o, kRepeatRow, kSectionRepeat, "repeat");

    std::vector<const RowResult*> rows;
    json rows_json = json::array();
    for (const RowResult& r : checklist) {
        rows.push_back(&r);
        rows_json.push_back(rowJson(r));
    }
    std::vector<const RowResult*> rows_with_repeat(rows.begin(), rows.end() - 1);
    rows_with_repeat.push_back(&repeat);

    json record;
    record["calibration_version"] = kCalibrationVersion;
    record["calibration_id"] = o.calibration_id;
    record["project_version"] = riptide::kProjectVersion;
    record["git_describe"] = riptide::kGitDescribe;
    // The thread count is left out: it does not change a number (D44).
    record["settings"] = {{"seed", o.seed},
                          {"sigma_Ep_MeV", o.resolution.sigma_Ep},
                          {"sigma_theta_rad", o.resolution.sigma_theta},
                          {"pilot", o.pilot},
                          {"levels", st::kDirectionCredibleLevels},
                          {"En_domain", domainJson(setup.domain)},
                          {"theta_points", st::kTrackAngleGridPoints},
                          {"theta_obs_points", st::kTrackObservedGridPoints},
                          {"sphere_pixels", st::kSphereGridPixels},
                          {"cap_pixels", st::kCapGridPixels},
                          {"cap_delta_log", st::kCapDeltaLog}};
    record["checklist"] = {{"rows", rows_json}, {"repeat", rowJson(repeat)}};
    record["contraction"] = {{"slope_min_events", kSlopeMinEvents},
                             {"slope", slopeOf(rows)},
                             {"slope_with_repeat", slopeOf(rows_with_repeat)}};
    record["prior"] = priorRows(setup, o, checklist);

    fs::create_directories(out_dir);
    riptide::io::writeFileAtomically(out_dir / kCalibrationName, record.dump(2) + "\n");
    std::cout << "riptide_case_b_calibrate: output in " << out_dir << '\n';
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    bool show_help = false;
    Options o;
    const auto cli =
        lyra::cli() | lyra::help(show_help) |
        lyra::opt(o.calibration_id, "id")["--calibration-id"]("name of the output directory")
            .required() |
        lyra::opt(o.data_dir, "dir")["--data-dir"]("data root (default: the environment)") |
        lyra::opt(o.seed, "n")["--seed"]("seed of the pseudo-experiments (default 20260907)") |
        lyra::opt(o.resolution.sigma_Ep, "MeV")["--sigma-Ep"]("Ep resolution (default 0.10)") |
        lyra::opt(o.resolution.sigma_theta,
                  "rad")["--sigma-theta"]("track resolution per component (default 0.08)") |
        lyra::opt(o.threads, "n")["--threads"]("threads across experiments (default 4)") |
        lyra::opt(o.pilot)["--pilot"]("counts / 10, same N: to measure time and memory");
    const auto result = cli.parse(lyra::args(argc, argv));
    if (show_help) {
        std::cout << cli << '\n';
        return kExitOk;
    }
    if (!result) {
        std::cerr << "riptide_case_b_calibrate: " << result.message() << "\n\n" << cli << '\n';
        return kExitUsage;
    }
    if (o.threads == 0) {
        std::cerr << "riptide_case_b_calibrate: --threads must be at least 1\n";
        return kExitUsage;
    }
    try {
        return calibrate(o);
    } catch (const std::exception& e) {
        std::cerr << "riptide_case_b_calibrate: " << e.what() << '\n';
        return kExitRuntime;
    }
}
