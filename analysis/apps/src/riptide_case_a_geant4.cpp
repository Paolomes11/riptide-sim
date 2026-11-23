// Case A on Geant4 events (D37, D38): reads one riptide_sim run, applies the
// fixed cut chain, converts the HP proton energies, draws the detector
// resolution with a recorded seed and splits the selected events into
// disjoint experiments, then runs the rows of the toy calibration on them.
// With --toy the pool comes from the toy generator instead, with the same
// truth and the same protocol, so the difference of the two outputs is the
// cost of the detector realism. With --acceptance the theta_p prior of the
// reference and wide-E_n models is the measured acceptance of a flat run
// (D40), and the toy draws its events with the same acceptance, so it stays
// the perfect-model twin of Geant4. It writes case_a.json under the stats
// directory of the data root; the file holds no wall time, so two runs with
// the same settings give the same bytes. Section timings go to stdout only.

#include <lyra/lyra.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
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
#include "riptide/apps/prepare_cli.hpp"
#include "riptide/common/data_paths.hpp"
#include "riptide/common/event_sample.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/random.hpp"
#include "riptide/common/version.hpp"
#include "riptide/io/acceptance_record.hpp"
#include "riptide/io/analysis_record.hpp"
#include "riptide/io/event_columns.hpp"
#include "riptide/io/response.hpp"
#include "riptide/io/root_reader.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"
#include "riptide/stats/acceptance.hpp"
#include "riptide/stats/calibration.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment.hpp"
#include "riptide/stats/split.hpp"
#include "riptide/stats/validate.hpp"

namespace {

namespace fs = std::filesystem;
namespace apps = riptide::apps;
namespace io = riptide::io;
namespace st = riptide::stats;
using Clock = std::chrono::steady_clock;
using nlohmann::json;

constexpr const char* kAppName = "riptide_case_a_geant4";
constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
constexpr std::uint64_t kProjectSeed = 20260907;
constexpr st::CaseAResolution kDefaultResolution{0.10, 0.08}; // MeV, rad
constexpr int kCaseAVersion = 2;
constexpr const char* kCaseAName = "case_a.json";
// Generous upper bound on the output: one small JSON file.
constexpr std::uintmax_t kOutputBytes = std::uintmax_t{1} << 20U; // 1 MiB

// The cut chain of Case A (D30, D34): one n-p scattering, a proton, and the
// proton contained, so its deposit is its energy.
constexpr const char* kCaseACuts = "single_np,has_proton,contained";

// The protocol of the toy calibration, so the rows can be compared one by one.
const std::vector<double> kLevels{0.68, 0.90, 0.95};
constexpr std::size_t kChecklistEvents = 20000;
constexpr std::size_t kChecklistBiasBins = 5;
constexpr double kWideEnMax = 30.0; // MeV: the alternative E_n prior on [10, 30]

struct SizedRow {
    std::size_t events;      // N per experiment
    std::size_t experiments; // planned M
};
const std::vector<SizedRow> kContractionRows{{10, 400}, {100, 200}, {1000, 40}};
const std::vector<SizedRow> kEnergyPriorRows{{1, 400}, {3, 400}, {10, 400}, {30, 200}, {100, 200}};
constexpr SizedRow kThetaPriorRow{1000, 20};

// Disjoint sections of the split engine, one per row, so no two rows share a
// permutation; within a row no event is used twice (D37).
constexpr std::uint64_t kSectionChecklist = 0;
constexpr std::uint64_t kSectionThetaPrior = 2;
constexpr std::uint64_t kSectionContraction = 16;
constexpr std::uint64_t kSectionEnergyPrior = 32;

struct Options {
    std::string data_dir;
    std::string stats_id;
    std::string run_id;
    std::string hp_check_id;
    bool toy = false;
    std::string truth;
    std::size_t events = 0;
    std::string acceptance_id;
    std::uint64_t seed = kProjectSeed;
    std::uint64_t split_seed = kProjectSeed;
    st::CaseAResolution resolution = kDefaultResolution;
    bool pilot = false;
};

// Exactly one source: a Geant4 run, or the toy with its truth and size.
void requireOneSource(const Options& o) {
    if (o.toy) {
        if (!o.run_id.empty() || !o.hp_check_id.empty()) {
            throw apps::CliError("--toy excludes --run-id and --hp-from-check");
        }
        if (o.truth.empty() || o.events == 0) {
            throw apps::CliError("--toy needs --truth and a positive --events");
        }
        return;
    }
    if (o.run_id.empty()) {
        throw apps::CliError("give either --run-id or --toy");
    }
    if (!o.truth.empty() || o.events != 0) {
        throw apps::CliError("--truth and --events go with --toy only");
    }
}

json domainJson(const st::EnergyDomain& d) {
    return {{"En_min_MeV", d.En_min}, {"En_max_MeV", d.En_max}, {"points", d.points}};
}

json truthJson(const io::SourceSpectrum& t) {
    json out = {{"mode", io::spectrumModeName(t.mode)}};
    if (t.mode == io::SpectrumMode::Mono) {
        out["En_MeV"] = t.En_MeV;
    } else {
        out["En_min_MeV"] = t.En_min_MeV;
        out["En_max_MeV"] = t.En_max_MeV;
    }
    return out;
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

// The theta_p acceptance read from acceptance.json, and where it came from.
struct Acceptance {
    std::optional<st::AcceptanceRatio> ratio;
    json info; // null without --acceptance
};

Acceptance loadAcceptance(const fs::path& root, const std::string& id) {
    Acceptance out;
    if (id.empty()) {
        return out;
    }
    const json record = io::readJsonFile(riptide::data::statsDir(root, id) /
                                         std::string{io::kAcceptanceRecordName});
    const io::AcceptanceTable table = io::acceptanceTableFromRecord(record);
    out.ratio.emplace(st::AcceptanceHistogram{table.En_edges, table.theta_edges, table.counts});
    out.info = {{"stats_id", id},
                {"acceptance_version", record.at("acceptance_version")},
                {"run_id", record.at("source").at("run_id")},
                {"half", record.at("half")},
                {"entries", record.at("entries")},
                {"En_bins", table.En_edges.size() - 1},
                {"theta_bins", table.theta_edges.size() - 1},
                {"ratio_max", out.ratio->maximum()}};
    return out;
}

// The theta_p prior of a model on En: physical by default, the measured
// acceptance with one row per E_n point otherwise (D40).
st::CaseAModel referenceModel(const st::Grid& En, const st::Grid& theta,
                              const st::AcceptanceRatio* acceptance, const st::CaseAResolution& r) {
    if (acceptance == nullptr) {
        return {En, theta, st::physicalRecoilAngleLogPrior(theta), r};
    }
    return {En, theta, st::acceptanceRecoilAngleLogPrior(*acceptance, En, theta), r};
}

// The models of the calibration: the reference one (default E_n grid,
// physical or acceptance theta_p prior) and the two alternatives of the
// robustness checks. The flat theta_p prior stays flat with --acceptance: it
// is the alternative the theta_p prior row compares against.
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

    Models(const st::CaseAResolution& r, const st::AcceptanceRatio* acceptance)
        : physical(referenceModel(En, theta, acceptance, r)),
          flat_theta(En, theta, st::flatRecoilAngleLogPrior(theta), r),
          wide_En(referenceModel(En_wide, theta, acceptance, r)) {}
};

// The pool of events and what it came from.
struct Pool {
    st::PseudoSample sample;
    io::SourceSpectrum truth;
    json source;
};

Pool geant4Pool(const fs::path& root, const Options& o) {
    const std::optional<std::string> check =
        o.hp_check_id.empty() ? std::nullopt : std::optional<std::string>{o.hp_check_id};
    const io::RunInfo run = io::readRunInfo(root, o.run_id, check);
    const io::Selection selection = io::Selection::parse(kCaseACuts);
    io::CutFlow flow(selection.cuts().size());
    Pool pool;
    pool.truth =
        io::sourceSpectrumFromRunRecord(io::readJsonFile(run.run_dir / io::kRunRecordName));
    {
        riptide::EventSample events;
        {
            // The full columns are needed only for the selection.
            const io::EventColumns columns = io::readEventColumns(run, 0);
            events = io::selectEvents(columns, selection, flow);
        }
        io::convertHpEnergies(events, run.hp.target_mass_MeV);
        io::applyResolution(events, io::Resolution{o.resolution.sigma_Ep, o.resolution.sigma_theta},
                            o.seed);
        pool.sample = apps::pseudoSampleFromEvents(events);
    }
    pool.source = {{"kind", "geant4"},
                   {"run_id", run.run_id},
                   {"run_seed", run.seed},
                   {"run_record_version", run.run_record_version},
                   {"schema_version", run.schema_version},
                   {"run_git_describe", run.git_describe},
                   {"hp_mass",
                    {{"source", io::hpMassSourceName(run.hp.source)},
                     {"check_id", run.hp.check_id},
                     {"target_mass_MeV", run.hp.target_mass_MeV}}},
                   {"cuts", selection.toString()},
                   {"events_read", flow.total()},
                   {"events_selected", flow.selected()},
                   {"cut_flow", io::cutsJson(selection, flow)}};
    return pool;
}

Pool toyPool(const Models& m, const st::AcceptanceRatio* acceptance, const Options& o) {
    Pool pool;
    pool.truth = apps::parseTruth(o.truth);
    if (pool.truth.mode == io::SpectrumMode::Mono) {
        pool.sample =
            acceptance == nullptr
                ? st::pseudoSampleAtEnergy(pool.truth.En_MeV, o.events, o.resolution, o.seed)
                : st::pseudoSampleAtEnergy(pool.truth.En_MeV, o.events, o.resolution, o.seed,
                                           *acceptance);
    } else {
        const st::EnergyDomain truth{pool.truth.En_min_MeV, pool.truth.En_max_MeV, m.domain.points};
        pool.sample =
            acceptance == nullptr
                ? st::pseudoSampleFromPrior(truth, o.events, o.resolution, o.seed)
                : st::pseudoSampleFromPrior(truth, o.events, o.resolution, o.seed, *acceptance);
    }
    pool.source = {{"kind", "toy"}, {"events", o.events}};
    return pool;
}

// The disjoint experiments of one row, or nothing when fewer than two fit:
// the row is then recorded as skipped, with what was planned and available.
struct RowSplit {
    std::size_t planned = 0;
    std::vector<std::vector<std::size_t>> experiments;
};

RowSplit splitRow(const st::PseudoSample& pool, SizedRow row, std::uint64_t section,
                  const Options& o) {
    RowSplit out;
    out.planned = apps::pilotScaled(row.experiments, o.pilot);
    const std::size_t m = apps::rowExperiments(row.experiments, pool.size(), row.events, o.pilot);
    if (m > 0) {
        riptide::random::Engine engine = riptide::random::perEventEngine(o.split_seed, section);
        out.experiments = st::disjointExperiments(pool.size(), row.events, m, engine);
    }
    return out;
}

json rowHeader(const st::PseudoSample& pool, SizedRow row, const RowSplit& split) {
    json out = {{"events", row.events},
                {"experiments_planned", split.planned},
                {"experiments", split.experiments.size()}};
    if (split.experiments.empty()) {
        out["experiments_available"] = st::maxDisjointExperiments(pool.size(), row.events);
        out["skipped"] = "too few events";
    } else {
        out["skipped"] = nullptr;
    }
    return out;
}

// Flat source: the single-event checklist on a random subset of the pool, with
// the bias per E_n bin.
json checklist(const Models& m, const Pool& pool, const Options& o) {
    const Stopwatch watch("checklist");
    const std::size_t n =
        std::min(apps::pilotScaled(kChecklistEvents, o.pilot), pool.sample.size());
    json out = {{"events_planned", apps::pilotScaled(kChecklistEvents, o.pilot)}, {"events", n}};
    if (n < apps::kMinExperiments) {
        out["skipped"] = "too few events";
        return out;
    }
    riptide::random::Engine engine =
        riptide::random::perEventEngine(o.split_seed, kSectionChecklist);
    const std::vector<std::vector<std::size_t>> subset =
        st::disjointExperiments(pool.sample.size(), n, 1, engine);
    const st::PseudoSample sample = st::pseudoSampleSubset(pool.sample, subset.front());
    st::PosteriorSummaries s = st::makeSummaries(kLevels);
    st::appendSingleEventSummaries(m.physical, m.log_prior_En, sample, s);
    out["skipped"] = nullptr;
    out["ensemble"] = ensembleJson(st::summariseEnsemble(s));
    out["bias_bins"] = biasJson(st::biasCurve(s.truth, s.mean, kChecklistBiasBins));
    watch.report();
    return out;
}

json contraction(const Models& m, const Pool& pool, const Options& o) {
    json rows = json::array();
    for (std::size_t r = 0; r < kContractionRows.size(); ++r) {
        const SizedRow row = kContractionRows[r];
        const Stopwatch watch("contraction N = " + std::to_string(row.events));
        const RowSplit split = splitRow(pool.sample, row, kSectionContraction + r, o);
        json out = rowHeader(pool.sample, row, split);
        if (!split.experiments.empty()) {
            st::PosteriorSummaries s = st::makeSummaries(kLevels);
            for (const std::vector<std::size_t>& indices : split.experiments) {
                st::appendCombinedSummary(m.physical, m.log_prior_En,
                                          st::pseudoSampleSubset(pool.sample, indices),
                                          pool.truth.En_MeV, s);
            }
            const st::EnsembleSummary ensemble = st::summariseEnsemble(s);
            const double root_n = std::sqrt(static_cast<double>(row.events));
            out["ensemble"] = ensembleJson(ensemble);
            out["std_sqrt_N"] = ensemble.mean_std * root_n;
            out["std_sqrt_N_error"] = ensemble.std_error * root_n;
        }
        rows.push_back(out);
        watch.report();
    }
    return {{"rows", rows}};
}

json thetaPrior(const Models& m, const Pool& pool, const Options& o) {
    const Stopwatch watch("theta_p prior");
    const RowSplit split = splitRow(pool.sample, kThetaPriorRow, kSectionThetaPrior, o);
    json out = rowHeader(pool.sample, kThetaPriorRow, split);
    if (!split.experiments.empty()) {
        st::PosteriorSummaries physical = st::makeSummaries(kLevels);
        st::PosteriorSummaries flat = st::makeSummaries(kLevels);
        for (const std::vector<std::size_t>& indices : split.experiments) {
            // Both priors see the very same events.
            const st::PseudoSample sample = st::pseudoSampleSubset(pool.sample, indices);
            st::appendCombinedSummary(m.physical, m.log_prior_En, sample, pool.truth.En_MeV,
                                      physical);
            st::appendCombinedSummary(m.flat_theta, m.log_prior_En, sample, pool.truth.En_MeV,
                                      flat);
        }
        out["physical"] = ensembleJson(st::summariseEnsemble(physical));
        out["flat"] = ensembleJson(st::summariseEnsemble(flat));
    }
    watch.report();
    return out;
}

json energyPrior(const Models& m, const Pool& pool, const Options& o) {
    json rows = json::array();
    for (std::size_t r = 0; r < kEnergyPriorRows.size(); ++r) {
        const SizedRow row = kEnergyPriorRows[r];
        const Stopwatch watch("E_n prior N = " + std::to_string(row.events));
        const RowSplit split = splitRow(pool.sample, row, kSectionEnergyPrior + r, o);
        json out = rowHeader(pool.sample, row, split);
        if (!split.experiments.empty()) {
            st::PosteriorSummaries reference = st::makeSummaries(kLevels);
            st::PosteriorSummaries alternative = st::makeSummaries(kLevels);
            for (const std::vector<std::size_t>& indices : split.experiments) {
                const st::PseudoSample sample = st::pseudoSampleSubset(pool.sample, indices);
                st::appendCombinedSummary(m.physical, m.log_prior_En, sample, pool.truth.En_MeV,
                                          reference);
                st::appendCombinedSummary(m.wide_En, m.log_prior_En_wide, sample, pool.truth.En_MeV,
                                          alternative);
            }
            const st::PriorShift shift = st::priorShift(reference, alternative);
            out["shift_mean"] = shift.shift_mean;
            out["shift_error"] = shift.shift_error;
            out["shift_max_abs"] = shift.shift_max_abs;
            out["std_ratio_mean"] = shift.std_ratio_mean;
            out["reference"] = ensembleJson(st::summariseEnsemble(reference));
            out["alternative"] = ensembleJson(st::summariseEnsemble(alternative));
        }
        rows.push_back(out);
        watch.report();
    }
    return {{"alternative_domain", domainJson(m.wide)}, {"rows", rows}};
}

int run(const Options& o) {
    namespace data = riptide::data;
    const std::optional<std::string> cli_root =
        o.data_dir.empty() ? std::nullopt : std::optional<std::string>{o.data_dir};
    const fs::path root = data::resolveDataRoot(cli_root);
    const fs::path out_dir = data::statsDir(root, o.stats_id);
    if (fs::exists(out_dir)) {
        throw data::DataPathError("stats directory already exists: " + out_dir.string());
    }
    data::requireWritableRoot(root, kOutputBytes, data::kFreeSpaceMargin, false);

    std::cout << kAppName << ' ' << riptide::kProjectVersion << " (" << riptide::kGitDescribe << ")"
              << (o.toy ? ", toy" : ", geant4") << (o.pilot ? ", pilot" : "") << '\n';
    const Acceptance acceptance = loadAcceptance(root, o.acceptance_id);
    const st::AcceptanceRatio* ratio = acceptance.ratio ? &*acceptance.ratio : nullptr;
    const Models models(o.resolution, ratio);
    const Pool pool = o.toy ? toyPool(models, ratio, o) : geant4Pool(root, o);
    std::cout << "  pool: " << pool.sample.size() << " events\n";

    json record;
    record["case_a_version"] = kCaseAVersion;
    record["stats_id"] = o.stats_id;
    record["project_version"] = riptide::kProjectVersion;
    record["git_describe"] = riptide::kGitDescribe;
    record["source"] = pool.source;
    record["truth"] = truthJson(pool.truth);
    record["settings"] = {{"seed", o.seed},
                          {"split_seed", o.split_seed},
                          {"sigma_Ep_MeV", o.resolution.sigma_Ep},
                          {"sigma_theta_rad", o.resolution.sigma_theta},
                          {"pilot", o.pilot},
                          {"levels", kLevels},
                          {"En_domain", domainJson(models.domain)},
                          {"theta_points", models.theta.size()},
                          {"acceptance", acceptance.info}};
    if (pool.truth.mode == io::SpectrumMode::Flat) {
        record["checklist"] = checklist(models, pool, o);
    } else {
        // A combined posterior assumes one E_n shared by every event.
        apps::requireMonoTruth(pool.sample, pool.truth.En_MeV);
        record["contraction"] = contraction(models, pool, o);
        record["theta_prior"] = thetaPrior(models, pool, o);
        record["energy_prior"] = energyPrior(models, pool, o);
    }

    fs::create_directories(out_dir);
    io::writeFileAtomically(out_dir / kCaseAName, record.dump(2) + "\n");
    std::cout << kAppName << ": output in " << out_dir << '\n';
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    bool show_help = false;
    Options o;
    const auto cli =
        lyra::cli() | lyra::help(show_help) |
        lyra::opt(o.stats_id, "id")["--stats-id"]("name of the output directory").required() |
        lyra::opt(o.data_dir, "dir")["--data-dir"]("data root (default: the environment)") |
        lyra::opt(o.run_id, "id")["--run-id"]("Geant4 source: the riptide_sim run") |
        lyra::opt(o.hp_check_id,
                  "id")["--hp-from-check"]("check with the HP mass (needed by run.json v1)") |
        lyra::opt(o.toy)["--toy"]("toy source (D38), with --truth and --events") |
        lyra::opt(o.truth, "spec")["--truth"]("mono:<En> or flat:<min>:<max>, MeV") |
        lyra::opt(o.events, "n")["--events"]("events of the toy pool") |
        lyra::opt(o.seed, "n")["--seed"]("seed of the smearing or the toy (default 20260907)") |
        lyra::opt(o.split_seed, "n")["--split-seed"]("seed of the split (default 20260907)") |
        lyra::opt(o.resolution.sigma_Ep, "MeV")["--sigma-Ep"]("Ep resolution (default 0.10)") |
        lyra::opt(o.resolution.sigma_theta,
                  "rad")["--sigma-theta"]("theta_p resolution (default 0.08)") |
        lyra::opt(o.acceptance_id,
                  "id")["--acceptance"]("stats id of acceptance.json: theta_p prior (D40)") |
        lyra::opt(o.pilot)["--pilot"]("counts / 10, same N: to measure time and memory");
    const auto result = cli.parse(lyra::args(argc, argv));
    if (show_help) {
        std::cout << cli << '\n';
        return kExitOk;
    }
    try {
        if (!result) {
            throw apps::CliError(result.message());
        }
        requireOneSource(o);
        if (o.toy) {
            static_cast<void>(apps::parseTruth(o.truth)); // a usage error, not a runtime one
        }
    } catch (const apps::CliError& e) {
        std::cerr << kAppName << ": " << e.what() << "\n\n" << cli << '\n';
        return kExitUsage;
    }
    try {
        return run(o);
    } catch (const std::exception& e) {
        std::cerr << kAppName << ": " << e.what() << '\n';
        return kExitRuntime;
    }
}
