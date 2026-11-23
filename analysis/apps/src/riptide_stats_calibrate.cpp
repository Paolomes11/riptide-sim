// Calibrates the Case A engine with pseudo-experiments drawn from the toy
// model (D35), so the port can be compared number by number with the toy
// before any Geant4 data: the checklist on single events, the contraction of
// the combined posterior with N and the robustness to the theta_p and E_n
// priors. It writes calibration.json under the stats directory of the data
// root; the file holds no wall time, so two runs with the same settings give
// the same bytes. Section timings go to stdout only.

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
#include "riptide/common/git_info.hpp"
#include "riptide/common/version.hpp"
#include "riptide/io/analysis_record.hpp"
#include "riptide/stats/calibration.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment.hpp"
#include "riptide/stats/validate.hpp"

namespace {

namespace fs = std::filesystem;
namespace st = riptide::stats;
using Clock = std::chrono::steady_clock;
using nlohmann::json;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
constexpr std::uint64_t kProjectSeed = 20260907;
constexpr st::CaseAResolution kToyResolution{0.10, 0.08}; // MeV, rad
constexpr int kCalibrationVersion = 1;
constexpr const char* kCalibrationName = "calibration.json";
// Generous upper bound on the output: one small JSON file.
constexpr std::uintmax_t kOutputBytes = std::uintmax_t{1} << 20U; // 1 MiB

// The protocol of the toy calibration, which the numbers are compared with.
const std::vector<double> kLevels{0.68, 0.90, 0.95};
constexpr std::size_t kChecklistEvents = 20000;
constexpr std::size_t kChecklistBiasBins = 5;
constexpr double kNarrowEnMin = 11.0;  // MeV: truth away from the prior edges
constexpr double kNarrowEnMax = 19.0;  // MeV
constexpr double kThetaPriorEn = 15.0; // MeV
constexpr std::size_t kThetaPriorEvents = 1000;
constexpr std::size_t kThetaPriorRepeats = 20;
constexpr double kWideEnMax = 30.0; // MeV: the alternative E_n prior on [10, 30]

struct SizedRow {
    std::size_t events;      // N per experiment
    std::size_t experiments; // M
};
const std::vector<SizedRow> kContractionRows{{10, 400}, {100, 200}, {1000, 40}};
const std::vector<SizedRow> kEnergyPriorRows{{1, 400}, {3, 400}, {10, 400}, {30, 200}, {100, 200}};

// Disjoint event-id sections, so no two draws of the calibration share an
// engine; a section of sized rows uses one section per row.
constexpr std::uint64_t kSectionChecklistPrior = 0;
constexpr std::uint64_t kSectionChecklistNarrow = 1;
constexpr std::uint64_t kSectionThetaPrior = 2;
constexpr std::uint64_t kSectionContraction = 16;
constexpr std::uint64_t kSectionEnergyPrior = 32;

struct Options {
    std::string data_dir;
    std::string calibration_id;
    std::uint64_t seed = kProjectSeed;
    st::CaseAResolution resolution = kToyResolution;
    bool pilot = false;
};

// The pilot divides every count by ten (at least two entries) but keeps N: it
// measures time and memory, not the physics.
std::size_t scaled(std::size_t count, bool pilot) {
    return riptide::apps::pilotScaled(count, pilot);
}

json domainJson(const st::EnergyDomain& d) {
    return {{"En_min_MeV", d.En_min}, {"En_max_MeV", d.En_max}, {"points", d.points}};
}

json ensembleJson(const st::EnsembleSummary& e) {
    return {{"entries", e.entries},     {"bias_mean", e.bias_mean},   {"bias_error", e.bias_error},
            {"rms_error", e.rms_error}, {"mean_std", e.mean_std},     {"std_error", e.std_error},
            {"pull_mean", e.pull_mean}, {"pull_width", e.pull_width}, {"levels", e.levels},
            {"coverage", e.coverage}};
}

json biasJson(const st::BiasCurve& b) {
    // An empty bin has a NaN bias, which JSON writes as null.
    return {{"centers", b.centers}, {"bias", b.bias}, {"counts", b.counts}};
}

class Stopwatch {
public:
    explicit Stopwatch(std::string name) : name_(std::move(name)) {}
    void report() const {
        const std::chrono::duration<double> elapsed = Clock::now() - start_;
        std::cout << "  " << name_ << ": " << elapsed.count() << " s\n" << std::flush;
    }

private:
    std::string name_;
    Clock::time_point start_ = Clock::now();
};

// The models of the calibration: the reference one (default E_n grid,
// physical theta_p prior) and the two alternatives of the robustness checks.
struct Models {
    st::EnergyDomain domain;
    st::EnergyDomain wide = st::sameStepDomain(domain, kWideEnMax);
    st::Grid En = st::energyGrid(domain);
    st::Grid En_wide = st::energyGrid(wide);
    st::Grid theta = st::recoilAngleGrid();
    std::vector<double> log_prior_En = st::flatEnergyLogPrior(En);
    std::vector<double> log_prior_En_wide = st::flatEnergyLogPrior(En_wide);
    st::CaseAModel physical;
    st::CaseAModel flat_theta;
    st::CaseAModel wide_En;

    explicit Models(const st::CaseAResolution& r)
        : physical(En, theta, st::physicalRecoilAngleLogPrior(theta), r),
          flat_theta(En, theta, st::flatRecoilAngleLogPrior(theta), r),
          wide_En(En_wide, theta, st::physicalRecoilAngleLogPrior(theta), r) {}
};

json checklist(const Models& m, const Options& o, const st::EnergyDomain& truth,
               std::uint64_t section, const std::string& name) {
    const Stopwatch watch("checklist " + name);
    const std::size_t n = scaled(kChecklistEvents, o.pilot);
    const st::PseudoSample sample = st::pseudoSampleFromPrior(
        truth, n, o.resolution, o.seed, st::experimentFirstEventId(section, 0));
    st::PosteriorSummaries s = st::makeSummaries(kLevels);
    st::appendSingleEventSummaries(m.physical, m.log_prior_En, sample, s);
    json out = {{"name", name},
                {"truth_domain", domainJson(truth)},
                {"events", n},
                {"ensemble", ensembleJson(st::summariseEnsemble(s))},
                {"bias_bins", biasJson(st::biasCurve(s.truth, s.mean, kChecklistBiasBins))}};
    watch.report();
    return out;
}

json contraction(const Models& m, const Options& o, const st::EnergyDomain& truth) {
    json rows = json::array();
    for (std::size_t r = 0; r < kContractionRows.size(); ++r) {
        const SizedRow row = kContractionRows[r];
        const Stopwatch watch("contraction N = " + std::to_string(row.events));
        const std::size_t experiments = scaled(row.experiments, o.pilot);
        st::PosteriorSummaries s = st::makeSummaries(kLevels);
        for (std::size_t k = 0; k < experiments; ++k) {
            const st::PseudoExperiment e = st::pseudoExperimentFromDomain(
                truth, row.events, o.resolution, o.seed, kSectionContraction + r, k);
            st::appendCombinedSummary(m.physical, m.log_prior_En, e.sample, e.truth_En, s);
        }
        const st::EnsembleSummary ensemble = st::summariseEnsemble(s);
        const double root_n = std::sqrt(static_cast<double>(row.events));
        rows.push_back({{"events", row.events},
                        {"experiments", experiments},
                        {"ensemble", ensembleJson(ensemble)},
                        {"std_sqrt_N", ensemble.mean_std * root_n},
                        {"std_sqrt_N_error", ensemble.std_error * root_n}});
        watch.report();
    }
    return {{"truth_domain", domainJson(truth)}, {"rows", rows}};
}

json thetaPrior(const Models& m, const Options& o) {
    const Stopwatch watch("theta_p prior");
    const std::size_t repeats = scaled(kThetaPriorRepeats, o.pilot);
    st::PosteriorSummaries physical = st::makeSummaries(kLevels);
    st::PosteriorSummaries flat = st::makeSummaries(kLevels);
    for (std::size_t k = 0; k < repeats; ++k) {
        // Both priors see the very same events.
        const st::PseudoExperiment e = st::pseudoExperimentAtEnergy(
            kThetaPriorEn, kThetaPriorEvents, o.resolution, o.seed, kSectionThetaPrior, k);
        st::appendCombinedSummary(m.physical, m.log_prior_En, e.sample, e.truth_En, physical);
        st::appendCombinedSummary(m.flat_theta, m.log_prior_En, e.sample, e.truth_En, flat);
    }
    json out = {{"En_MeV", kThetaPriorEn},
                {"events", kThetaPriorEvents},
                {"experiments", repeats},
                {"physical", ensembleJson(st::summariseEnsemble(physical))},
                {"flat", ensembleJson(st::summariseEnsemble(flat))}};
    watch.report();
    return out;
}

json energyPrior(const Models& m, const Options& o) {
    json rows = json::array();
    for (std::size_t r = 0; r < kEnergyPriorRows.size(); ++r) {
        const SizedRow row = kEnergyPriorRows[r];
        const Stopwatch watch("E_n prior N = " + std::to_string(row.events));
        const std::size_t experiments = scaled(row.experiments, o.pilot);
        st::PosteriorSummaries reference = st::makeSummaries(kLevels);
        st::PosteriorSummaries alternative = st::makeSummaries(kLevels);
        for (std::size_t k = 0; k < experiments; ++k) {
            const st::PseudoExperiment e = st::pseudoExperimentFromDomain(
                m.domain, row.events, o.resolution, o.seed, kSectionEnergyPrior + r, k);
            st::appendCombinedSummary(m.physical, m.log_prior_En, e.sample, e.truth_En, reference);
            st::appendCombinedSummary(m.wide_En, m.log_prior_En_wide, e.sample, e.truth_En,
                                      alternative);
        }
        const st::PriorShift shift = st::priorShift(reference, alternative);
        rows.push_back({{"events", row.events},
                        {"experiments", experiments},
                        {"shift_mean", shift.shift_mean},
                        {"shift_error", shift.shift_error},
                        {"shift_max_abs", shift.shift_max_abs},
                        {"std_ratio_mean", shift.std_ratio_mean},
                        {"reference", ensembleJson(st::summariseEnsemble(reference))},
                        {"alternative", ensembleJson(st::summariseEnsemble(alternative))}});
        watch.report();
    }
    return {{"truth_domain", domainJson(m.domain)},
            {"alternative_domain", domainJson(m.wide)},
            {"rows", rows}};
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

    std::cout << "riptide_stats_calibrate " << riptide::kProjectVersion << " ("
              << riptide::kGitDescribe << ")" << (o.pilot ? ", pilot" : "") << '\n';
    const Models models(o.resolution);
    const st::EnergyDomain narrow{kNarrowEnMin, kNarrowEnMax, models.domain.points};

    json record;
    record["calibration_version"] = kCalibrationVersion;
    record["calibration_id"] = o.calibration_id;
    record["project_version"] = riptide::kProjectVersion;
    record["git_describe"] = riptide::kGitDescribe;
    record["settings"] = {{"seed", o.seed},
                          {"sigma_Ep_MeV", o.resolution.sigma_Ep},
                          {"sigma_theta_rad", o.resolution.sigma_theta},
                          {"pilot", o.pilot},
                          {"levels", kLevels},
                          {"En_domain", domainJson(models.domain)},
                          {"theta_points", models.theta.size()}};
    record["checklist"] = json::array(
        {checklist(models, o, models.domain, kSectionChecklistPrior, "truth_from_prior"),
         checklist(models, o, narrow, kSectionChecklistNarrow, "truth_narrow")});
    record["contraction"] = contraction(models, o, narrow);
    record["theta_prior"] = thetaPrior(models, o);
    record["energy_prior"] = energyPrior(models, o);

    fs::create_directories(out_dir);
    riptide::io::writeFileAtomically(out_dir / kCalibrationName, record.dump(2) + "\n");
    std::cout << "riptide_stats_calibrate: output in " << out_dir << '\n';
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
                  "rad")["--sigma-theta"]("theta_p resolution (default 0.08)") |
        lyra::opt(o.pilot)["--pilot"]("counts / 10, same N: to measure time and memory");
    const auto result = cli.parse(lyra::args(argc, argv));
    if (show_help) {
        std::cout << cli << '\n';
        return kExitOk;
    }
    if (!result) {
        std::cerr << "riptide_stats_calibrate: " << result.message() << "\n\n" << cli << '\n';
        return kExitUsage;
    }
    try {
        return calibrate(o);
    } catch (const std::exception& e) {
        std::cerr << "riptide_stats_calibrate: " << e.what() << '\n';
        return kExitRuntime;
    }
}
