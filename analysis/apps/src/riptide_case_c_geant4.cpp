// Case C on Geant4 events (D48): reads one flat riptide_sim run with a fixed
// neutron direction, applies the cut chain of Case A and B, converts the HP
// proton energies and draws the detector resolution with a recorded seed.
// The pool has the flat source spectrum shaped by the cross section and the
// selection; each experiment keeps pool events with probability
// f / (c g_hat), so its E_n follow the truncated Gaussian of its truth. The
// rows are those of the toy checklist (D47) plus the repeat at N = 1000,
// with the truths (mu_E, sigma_E) of the 8a calibration at the direction of
// the run. Three modes, run in this order:
//   --pool-only   the resampling alone, no reconstruction: how many
//                 experiments each row gets from the pool;
//   --toy --match the toy twin of a pool record: same truths, N and M, E_n
//                 drawn from f, the kinematics and response of the toy;
//   (default)     the reconstruction of the Geant4 experiments, with --twin
//                 naming the toy record the rms values are compared with.
// Every number with a threshold of the written predictions (8b) carries its
// verdict next to it, and the verdicts are tallied. Experiments run on
// several threads; each has its own engines and output slot, so the numbers
// do not depend on --threads. It writes case_c.json under the stats
// directory of the data root; the file holds no wall time.

#include <lyra/lyra.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "riptide/apps/case_a_geant4.hpp"
#include "riptide/apps/case_b_geant4.hpp"
#include "riptide/apps/prepare_cli.hpp"
#include "riptide/common/data_paths.hpp"
#include "riptide/common/direction.hpp"
#include "riptide/common/event_sample.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"
#include "riptide/common/version.hpp"
#include "riptide/io/acceptance_record.hpp"
#include "riptide/io/analysis_record.hpp"
#include "riptide/io/event_columns.hpp"
#include "riptide/io/response.hpp"
#include "riptide/io/root_reader.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"
#include "riptide/stats/acceptance.hpp"
#include "riptide/stats/calibration_b.hpp"
#include "riptide/stats/calibration_c.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_c.hpp"
#include "riptide/stats/case_c_posterior.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"
#include "riptide/stats/pseudo_experiment_c.hpp"
#include "riptide/stats/resample_c.hpp"
#include "riptide/stats/split.hpp"

namespace {

namespace fs = std::filesystem;
namespace apps = riptide::apps;
namespace io = riptide::io;
namespace st = riptide::stats;
namespace units = riptide::units;
using riptide::geometry::Direction;
using Clock = std::chrono::steady_clock;
using nlohmann::json;

constexpr const char* kAppName = "riptide_case_c_geant4";
constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
constexpr std::uint64_t kProjectSeed = 20260907;
constexpr st::CaseBResolution kDefaultResolution{0.10, 0.08}; // MeV, rad
constexpr unsigned kDefaultThreads = 4;
constexpr int kCaseCGeant4Version = 1;
constexpr const char* kCaseCName = "case_c.json";
// Generous upper bound on the output: one small JSON file.
constexpr std::uintmax_t kOutputBytes = std::uintmax_t{1} << 20U; // 1 MiB

// The cut chain of Case A and B (D30, D34, D45).
constexpr const char* kCaseCCuts = "single_np,has_proton,contained";

// The flat runs of 7b reach 19.5 MeV only, so the truncation, the E_n grid
// and the mu_E grid all stop there (D48).
constexpr double kDomainMin = 10.0 * units::MeV;
constexpr double kDomainMax = 19.5 * units::MeV;
// The run spectrum must be the domain: the same bounds up to rounding.
constexpr double kSpectrumRelTolerance = 1e-12;

// The rows: the toy checklist of D47 on the sections of the 8a calibration,
// so the truths are those of 8a with omega_n replaced by the run direction;
// the last one is the repeat at N = 1000, on which the verdicts at N = 1000
// are read. Rows share the pool (D48): each walks its own permutation.
struct RowPlan {
    std::size_t events;      // N per experiment
    std::size_t experiments; // planned M
    std::uint64_t section;   // of the truths, the split and the twin events
    bool repeat;
};
const std::vector<RowPlan> kRows{{50, 200, 16, false},
                                 {150, 200, 17, false},
                                 {300, 100, 18, false},
                                 {1000, 40, 19, false},
                                 {1000, 240, 24, true}};

// Thresholds of the predictions, section 4.
constexpr double kBiasStdFraction = 0.5;     // |bias| < this x mean posterior std
constexpr std::size_t kBiasEvents = 1000;    // the bias is judged at this N
constexpr double kDirectionBiasSigmas = 2.0; // axis runs: each component within this
constexpr double kObliqueBiasFraction = 0.5; // oblique run: norm < this x sigma_c
constexpr double kRmsRatioMax = 1.15;        // Geant4 / twin
constexpr double kPullTolerance = 0.1;       // pull width and pull rms 1 +- this
constexpr double kCoverageSigmas = 3.0;      // binomial
// The declared exception: the log sigma_E coverage at N = 50 is low in the
// toy too; here it is compared with the twin, without a verdict that counts.
constexpr std::size_t kLowCoverageEvents = 50;
// An axis run has one component of +-1; the bias rule of the oblique run
// applies otherwise.
constexpr double kAxisTolerance = 1e-9;

enum class Mode { pool_only, toy, geant4 };

const char* modeName(Mode m) {
    switch (m) {
    case Mode::pool_only:
        return "pool_only";
    case Mode::toy:
        return "toy";
    case Mode::geant4:
        return "geant4";
    }
    return "geant4";
}

struct Options {
    std::string data_dir;
    std::string stats_id;
    std::string run_id;
    std::string hp_check_id;
    bool pool_only = false;
    bool toy = false;
    std::string match_id;
    std::string twin_id;
    std::string acceptance_id;
    std::uint64_t seed = kProjectSeed;
    std::uint64_t split_seed = kProjectSeed;
    st::CaseBResolution resolution = kDefaultResolution;
    unsigned threads = kDefaultThreads;
    bool pilot = false;

    [[nodiscard]] Mode mode() const noexcept {
        if (pool_only) {
            return Mode::pool_only;
        }
        return toy ? Mode::toy : Mode::geant4;
    }
};

void requireConsistent(const Options& o) {
    if (o.pool_only && o.toy) {
        throw apps::CliError("--pool-only and --toy exclude each other");
    }
    if (o.toy) {
        if (!o.run_id.empty() || !o.hp_check_id.empty() || !o.twin_id.empty()) {
            throw apps::CliError("--toy excludes --run-id, --hp-from-check and --twin");
        }
        if (o.match_id.empty()) {
            throw apps::CliError("--toy needs --match, the pool record it is the twin of");
        }
        return;
    }
    if (o.run_id.empty()) {
        throw apps::CliError("give --run-id, or --toy with --match");
    }
    if (!o.match_id.empty()) {
        throw apps::CliError("--match goes with --toy only");
    }
    if (o.pool_only) {
        if (!o.twin_id.empty() || !o.acceptance_id.empty()) {
            throw apps::CliError("--pool-only excludes --twin and --acceptance");
        }
        return;
    }
    if (o.twin_id.empty()) {
        throw apps::CliError("the Geant4 reconstruction needs --twin, its toy record");
    }
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

// The role of a verdict: `counted` items make the false-alarm statistics; a
// `reference` item is a comparison reported without a verdict that counts.
enum class Role { counted, reference };

const char* roleName(Role r) {
    return r == Role::counted ? "counted" : "reference";
}

struct Tally {
    std::size_t counted = 0;
    std::vector<std::string> counted_failures;
    std::vector<std::string> reference_failures;

    json check(const std::string& name, double value, double lo, double hi, Role role) {
        const bool pass = value >= lo && value <= hi;
        if (role == Role::counted) {
            ++counted;
        }
        if (!pass) {
            (role == Role::counted ? counted_failures : reference_failures).push_back(name);
        }
        return {{"value", value}, {"lo", lo}, {"hi", hi}, {"pass", pass}, {"role", roleName(role)}};
    }

    [[nodiscard]] json summary() const {
        return {{"counted", counted},
                {"counted_failures", counted_failures},
                {"reference_failures", reference_failures}};
    }
};

double binomialBound(double p, std::size_t n) {
    return kCoverageSigmas * std::sqrt(p * (1.0 - p) / static_cast<double>(n));
}

json directionJson(const Direction& d) {
    return json::array({d.x, d.y, d.z});
}

Direction directionFromJson(const json& j) {
    return {j.at(0).get<double>(), j.at(1).get<double>(), j.at(2).get<double>()};
}

bool isAxis(const Direction& d) {
    const double largest = std::max({std::abs(d.x), std::abs(d.y), std::abs(d.z)});
    return largest >= 1.0 - kAxisTolerance;
}

json domainJson(const st::EnergyDomain& d) {
    return {{"En_min_MeV", d.En_min}, {"En_max_MeV", d.En_max}, {"points", d.points}};
}

json hyperJson(const st::HyperEnsemble& e) {
    return {{"entries", e.entries},
            {"bias", e.bias},
            {"bias_error", e.bias_error},
            {"rms", e.rms},
            {"std_mean", e.std_mean},
            {"pull_mean", e.pull_mean},
            {"pull_mean_error", e.pull_mean_error},
            {"pull_width", e.pull_width},
            {"coverage", e.coverage}};
}

json omegaJson(const st::CaseBEnsemble& e, const st::DirectionBias& b) {
    return {{"entries", e.entries},
            {"error_mean_rad", e.error_mean},
            {"error_mean_error_rad", e.error_mean_error},
            {"error_rms_rad", e.error_rms},
            {"error_rms_deg", e.error_rms / units::deg},
            {"sigma_mean_rad", e.sigma_mean},
            {"sigma_mean_deg", e.sigma_mean / units::deg},
            {"pull_rms", e.pull_rms},
            {"scaled_pull_median", e.scaled_pull_median},
            {"scaled_pull_tail_fraction", e.scaled_pull_tail_fraction},
            {"truth_in_cap_fraction", e.truth_in_cap_fraction},
            {"coverage", e.coverage},
            {"tangent_bias_rad", b.mean},
            {"tangent_bias_error_rad", b.error},
            {"tangent_bias_norm_rad", std::hypot(b.mean[0], b.mean[1])}};
}

json ensembleJson(const st::CaseCEnsemble& e) {
    return {{"mu_MeV", hyperJson(e.mu)},
            {"log_sigma", hyperJson(e.log_sigma)},
            {"omega", omegaJson(e.direction, e.direction_bias)}};
}

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

st::EnergyDomain caseCDomain() {
    return {kDomainMin, kDomainMax, st::kEnergyGridPoints};
}

// The model of 8a on the domain of 8b; with --acceptance the selection
// weights of D40 on the same grids.
struct Setup {
    st::EnergyDomain domain = caseCDomain();
    st::Grid En = st::energyGrid(domain);
    st::Grid theta = st::recoilAngleGrid(st::kTrackAngleGridPoints);
    st::CaseCModel model;
    std::vector<Direction> coarse_grid = st::sphereGrid(st::kSphereGridPixels);

    Setup(const st::CaseBResolution& r, const st::AcceptanceRatio* acceptance)
        : model(En, theta, st::Grid{0.0, units::pi, st::kTrackObservedGridPoints}, r,
                acceptance == nullptr ? std::vector<double>{}
                                      : st::caseBAcceptanceLogWeights(*acceptance, En, theta)) {}
};

// The selected Geant4 events and their true direction.
struct Pool {
    st::CaseBPseudoSample sample;
    Direction omega_n;
    json source;
};

void requireDomainSpectrum(const io::SourceSpectrum& s) {
    const auto close = [](double a, double b) {
        return std::abs(a - b) <= kSpectrumRelTolerance * std::abs(b);
    };
    if (s.mode != io::SpectrumMode::Flat || !close(s.En_min_MeV, kDomainMin) ||
        !close(s.En_max_MeV, kDomainMax)) {
        throw std::invalid_argument("the run must have the flat spectrum on the Case C domain [" +
                                    std::to_string(kDomainMin) + ", " + std::to_string(kDomainMax) +
                                    "] MeV");
    }
}

Pool geant4Pool(const fs::path& root, const Options& o) {
    const std::optional<std::string> check =
        o.hp_check_id.empty() ? std::nullopt : std::optional<std::string>{o.hp_check_id};
    const io::RunInfo run = io::readRunInfo(root, o.run_id, check);
    const json run_record = io::readJsonFile(run.run_dir / io::kRunRecordName);
    requireDomainSpectrum(io::sourceSpectrumFromRunRecord(run_record));
    const io::Selection selection = io::Selection::parse(kCaseCCuts);
    io::CutFlow flow(selection.cuts().size());
    Pool pool;
    pool.omega_n = apps::fixedDirectionFromRunRecord(run_record);
    {
        riptide::EventSample events;
        {
            // The full columns are needed only for the selection.
            const io::EventColumns columns = io::readEventColumns(run, 0);
            events = io::selectEvents(columns, selection, flow);
        }
        apps::requireFixedDirection(events, pool.omega_n);
        io::convertHpEnergies(events, run.hp.target_mass_MeV);
        io::applyResolution(events, io::Resolution{o.resolution.sigma_Ep, o.resolution.sigma_theta},
                            o.seed);
        pool.sample = apps::caseBSampleFromEvents(events);
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

json poolDensityJson(const st::BinnedDensity& g) {
    const auto [lo, hi] = std::minmax_element(g.density.begin(), g.density.end());
    return {{"bin_width_MeV", g.step},
            {"bins", g.density.size()},
            {"density_min_per_MeV", *lo},
            {"density_max_per_MeV", *hi}};
}

std::vector<st::CaseCTruth> rowTruths(const RowPlan& p, std::size_t planned,
                                      const Direction& omega_n, const Options& o) {
    std::vector<st::CaseCTruth> truths;
    truths.reserve(planned);
    for (std::size_t k = 0; k < planned; ++k) {
        truths.push_back(
            st::caseCTruthAtDirection(st::kCaseCChecklistTruth, omega_n, o.seed, p.section, k));
    }
    return truths;
}

json resamplingJson(const st::ResampledRow& r, std::size_t events) {
    const std::size_t visited = std::accumulate(r.visited.begin(), r.visited.end(), std::size_t{0});
    json out = {{"exhausted", r.exhausted()}, {"visited_total", r.visited_total}};
    if (r.experiments.empty()) {
        out["visited_mean"] = nullptr;
        out["kept_fraction"] = nullptr;
        out["bound_mean"] = nullptr;
        out["bound_max"] = nullptr;
        return out;
    }
    const auto m = static_cast<double>(r.experiments.size());
    out["visited_mean"] = static_cast<double>(visited) / m;
    out["kept_fraction"] = m * static_cast<double>(events) / static_cast<double>(visited);
    out["bound_mean"] = std::accumulate(r.bounds.begin(), r.bounds.end(), 0.0) / m;
    out["bound_max"] = *std::max_element(r.bounds.begin(), r.bounds.end());
    return out;
}

// What the toy twin compares with, read back from its record.
struct TwinRow {
    std::size_t experiments = 0;
    double rms_mu = 0.0;
    double rms_log_sigma = 0.0;
    double rms_omega = 0.0;   // rad
    double sigma_omega = 0.0; // rad, mean declared resolution
    std::array<double, st::kHyperCredibleLevels.size()> coverage_log_sigma{};
};

std::string label(const RowPlan& p, std::size_t m, const std::string& what) {
    return std::string(p.repeat ? "repeat " : "") + "N=" + std::to_string(p.events) +
           " M=" + std::to_string(m) + " " + what;
}

const std::array<const char*, 3> kLevelNames{"0.68", "0.90", "0.95"};

// The section-4 items of one Geant4 row against its twin.
json geant4Checks(Tally& t, const RowPlan& p, const st::CaseCEnsemble& e, const TwinRow& twin,
                  const Direction& omega_n, Role role) {
    const std::size_t m = e.mu.entries;
    const auto name = [&](const std::string& what) { return label(p, m, what); };
    json out;
    const double pull_lo = 1.0 - kPullTolerance;
    const double pull_hi = 1.0 + kPullTolerance;
    out["pull_width_mu"] = t.check(name("pull width mu"), e.mu.pull_width, pull_lo, pull_hi, role);
    out["pull_width_log_sigma"] =
        t.check(name("pull width log sigma"), e.log_sigma.pull_width, pull_lo, pull_hi, role);
    out["pull_rms_omega"] =
        t.check(name("pull rms omega"), e.direction.pull_rms, pull_lo, pull_hi, role);
    const bool low_n = p.events == kLowCoverageEvents;
    json coverage = json::array();
    for (std::size_t l = 0; l < st::kHyperCredibleLevels.size(); ++l) {
        const double level = st::kHyperCredibleLevels[l];
        const double bound = binomialBound(level, m);
        const std::string tag = kLevelNames[l];
        json c;
        c["level"] = level;
        c["mu"] = t.check(name("coverage mu " + tag), e.mu.coverage[l], level - bound,
                          level + bound, role);
        if (low_n) {
            // Against the twin: both coverages carry a binomial error.
            const double ref = twin.coverage_log_sigma[l];
            const double ref_bound =
                kCoverageSigmas * std::sqrt(ref * (1.0 - ref) * (2.0 / static_cast<double>(m)));
            c["log_sigma"] =
                t.check(name("coverage log sigma vs twin " + tag), e.log_sigma.coverage[l],
                        ref - ref_bound, ref + ref_bound, Role::reference);
        } else {
            c["log_sigma"] = t.check(name("coverage log sigma " + tag), e.log_sigma.coverage[l],
                                     level - bound, level + bound, role);
        }
        c["omega"] = t.check(name("coverage omega " + tag), e.direction.coverage[l], level - bound,
                             level + bound, role);
        coverage.push_back(c);
    }
    out["coverage"] = coverage;
    out["rms_ratio_mu"] =
        t.check(name("rms ratio mu"), e.mu.rms / twin.rms_mu, 0.0, kRmsRatioMax, role);
    out["rms_ratio_log_sigma"] = t.check(
        name("rms ratio log sigma"), e.log_sigma.rms / twin.rms_log_sigma, 0.0, kRmsRatioMax, role);
    out["rms_ratio_omega"] = t.check(
        name("rms ratio omega"), e.direction.error_rms / twin.rms_omega, 0.0, kRmsRatioMax, role);
    // Reported without a threshold, as in 7b.
    out["sigma_ratio_omega"] = e.direction.sigma_mean / twin.sigma_omega;
    out["truth_in_cap"] =
        t.check(name("truth in cap"), e.direction.truth_in_cap_fraction, 1.0, 1.0, Role::reference);
    if (p.events != kBiasEvents) {
        return out;
    }
    const double mu_bound = kBiasStdFraction * e.mu.std_mean;
    const double ls_bound = kBiasStdFraction * e.log_sigma.std_mean;
    out["bias_mu"] = t.check(name("bias mu"), e.mu.bias, -mu_bound, mu_bound, role);
    out["bias_log_sigma"] =
        t.check(name("bias log sigma"), e.log_sigma.bias, -ls_bound, ls_bound, role);
    const st::DirectionBias& b = e.direction_bias;
    if (isAxis(omega_n)) {
        // On an axis the bias vanishes by symmetry: each tangent component
        // within kDirectionBiasSigmas of its error.
        json components = json::array();
        for (std::size_t c = 0; c < b.mean.size(); ++c) {
            const double bound = kDirectionBiasSigmas * b.error[c];
            components.push_back(t.check(name("bias omega component " + std::to_string(c)),
                                         b.mean[c], -bound, bound, role));
        }
        out["bias_omega"] = {{"rule", "each tangent component within 2 errors"},
                             {"components", components}};
    } else {
        const double bound = kObliqueBiasFraction * e.direction.sigma_mean;
        out["bias_omega"] = {{"rule", "norm of the mean tangent vector < 0.5 sigma_c"},
                             {"norm", t.check(name("bias omega norm"),
                                              std::hypot(b.mean[0], b.mean[1]), 0.0, bound, role)}};
    }
    return out;
}

// One row: the plan, how many experiments it has, and what came out.
struct RowRecord {
    std::size_t planned = 0;
    std::size_t experiments = 0;
    std::optional<st::CaseCEnsemble> ensemble; // empty if skipped or pool only
    json resampling;                           // null for the toy
};

json rowJson(const RowPlan& p, const RowRecord& r) {
    json out = {{"events", p.events},
                {"section", p.section},
                {"repeat", p.repeat},
                {"experiments_planned", r.planned},
                {"experiments", r.experiments},
                {"skipped", r.experiments < apps::kMinExperiments ? json("too few experiments")
                                                                  : json(nullptr)},
                {"resampling", r.resampling}};
    out["ensemble"] = r.ensemble ? ensembleJson(*r.ensemble) : json(nullptr);
    return out;
}

std::vector<st::CaseCExperimentResult> analyseGeant4Row(const Setup& s, const Pool& pool,
                                                        const std::vector<st::CaseCTruth>& truths,
                                                        const st::ResampledRow& row,
                                                        unsigned threads) {
    std::vector<st::CaseCExperimentResult> results(row.experiments.size());
    st::parallelForEach(results.size(), threads, [&](std::size_t k) {
        const st::CaseCPseudoExperiment experiment{
            truths[k], st::caseBPseudoSampleSubset(pool.sample, row.experiments[k])};
        results[k] = st::analyseCaseCPseudoExperiment(s.model, experiment, s.coarse_grid, false);
    });
    return results;
}

std::vector<st::CaseCExperimentResult>
analyseToyRow(const Setup& s, const RowPlan& p, const std::vector<st::CaseCTruth>& truths,
              std::size_t m, const st::AcceptanceRatio* acceptance, const Options& o) {
    std::vector<st::CaseCExperimentResult> results(m);
    st::parallelForEach(m, o.threads, [&](std::size_t k) {
        const st::CaseCPseudoExperiment experiment = st::caseCPseudoExperimentAtTruth(
            s.domain, p.events, o.resolution, truths[k], o.seed, p.section, k, acceptance);
        results[k] = st::analyseCaseCPseudoExperiment(s.model, experiment, s.coarse_grid, false);
    });
    return results;
}

// A record this one depends on must have used the same truths and split.
void requireSameSettings(const json& other, const Options& o, const std::string& id) {
    const json& s = other.at("settings");
    if (s.at("seed").get<std::uint64_t>() != o.seed ||
        s.at("split_seed").get<std::uint64_t>() != o.split_seed ||
        s.at("pilot").get<bool>() != o.pilot) {
        throw std::invalid_argument("record " + id +
                                    " has another --seed, --split-seed or --pilot");
    }
    if (other.at("rows").size() != kRows.size()) {
        throw std::invalid_argument("record " + id + " has another set of rows");
    }
    for (std::size_t r = 0; r < kRows.size(); ++r) {
        const json& row = other.at("rows").at(r);
        if (row.at("events").get<std::size_t>() != kRows[r].events ||
            row.at("section").get<std::uint64_t>() != kRows[r].section) {
            throw std::invalid_argument("record " + id + " has another row " + std::to_string(r));
        }
    }
}

json readCaseCRecord(const fs::path& root, const std::string& id, const char* kind) {
    json record = io::readJsonFile(riptide::data::statsDir(root, id) / kCaseCName);
    if (record.at("mode").get<std::string>() != kind) {
        throw std::invalid_argument("record " + id + " is not a " + std::string(kind) + " record");
    }
    return record;
}

TwinRow twinRow(const json& row) {
    const json& e = row.at("ensemble");
    TwinRow out;
    out.experiments = row.at("experiments").get<std::size_t>();
    out.rms_mu = e.at("mu_MeV").at("rms").get<double>();
    out.rms_log_sigma = e.at("log_sigma").at("rms").get<double>();
    out.rms_omega = e.at("omega").at("error_rms_rad").get<double>();
    out.sigma_omega = e.at("omega").at("sigma_mean_rad").get<double>();
    out.coverage_log_sigma =
        e.at("log_sigma").at("coverage").get<decltype(out.coverage_log_sigma)>();
    return out;
}

json recordHeader(const Options& o, const st::EnergyDomain& domain) {
    json record;
    record["case_c_geant4_version"] = kCaseCGeant4Version;
    record["stats_id"] = o.stats_id;
    record["mode"] = modeName(o.mode());
    record["project_version"] = riptide::kProjectVersion;
    record["git_describe"] = riptide::kGitDescribe;
    // The thread count is left out: it does not change a number.
    record["settings"] = {{"seed", o.seed},
                          {"split_seed", o.split_seed},
                          {"sigma_Ep_MeV", o.resolution.sigma_Ep},
                          {"sigma_theta_rad", o.resolution.sigma_theta},
                          {"pilot", o.pilot},
                          {"levels", st::kHyperCredibleLevels},
                          {"En_domain", domainJson(domain)},
                          {"truth_range",
                           {{"mu_min_MeV", st::kCaseCChecklistTruth.mu_min},
                            {"mu_max_MeV", st::kCaseCChecklistTruth.mu_max},
                            {"sigma_min_MeV", st::kCaseCChecklistTruth.sigma_min},
                            {"sigma_max_MeV", st::kCaseCChecklistTruth.sigma_max}}},
                          {"theta_points", st::kTrackAngleGridPoints},
                          {"theta_obs_points", st::kTrackObservedGridPoints},
                          {"sphere_pixels", st::kSphereGridPixels},
                          {"cap_pixels", st::kCapGridPixels},
                          {"cap_delta_log", st::kCapDeltaLog},
                          {"hyper_global_points", {st::kHyperMuPoints, st::kHyperSigmaPoints}},
                          {"hyper_sigma_range_MeV", {st::kHyperSigmaMin, st::kHyperSigmaMax}},
                          {"hyper_fine_points", st::kHyperFinePoints},
                          {"hyper_window_delta_log", st::kHyperWindowDeltaLog},
                          {"resample_bin_width_MeV", st::kResampleBinWidth}};
    return record;
}

std::string rowName(const RowPlan& p, std::size_t m) {
    return std::string(p.repeat ? "repeat" : "row") + " N = " + std::to_string(p.events) +
           ", M = " + std::to_string(m);
}

// Pool only, or the Geant4 reconstruction: the same resampling of each row.
json geant4Record(const fs::path& root, const Options& o) {
    const Pool pool = geant4Pool(root, o);
    std::cout << "  pool: " << pool.sample.size() << " events\n";
    const st::BinnedDensity g = st::histogramDensity(pool.sample.truth_En, caseCDomain());

    std::optional<json> twin;
    if (o.mode() == Mode::geant4) {
        twin = readCaseCRecord(root, o.twin_id, modeName(Mode::toy));
        requireSameSettings(*twin, o, o.twin_id);
        if (twin->at("source").at("match").at("run_id").get<std::string>() != o.run_id) {
            throw std::invalid_argument("the twin " + o.twin_id + " matches another run");
        }
        // The twin must share the model: an ideal Geant4 run against an
        // ideal toy, a weighted one against a toy drawn and fitted with the
        // same acceptance.
        const json& twin_acceptance = twin->at("settings").at("acceptance");
        const std::string twin_acceptance_id =
            twin_acceptance.is_null() ? std::string{}
                                      : twin_acceptance.at("stats_id").get<std::string>();
        if (twin_acceptance_id != o.acceptance_id) {
            throw std::invalid_argument("the twin " + o.twin_id + " has another --acceptance");
        }
    }
    const Acceptance acceptance = loadAcceptance(root, o.acceptance_id);
    std::optional<Setup> setup;
    if (o.mode() == Mode::geant4) {
        setup.emplace(o.resolution, acceptance.ratio ? &*acceptance.ratio : nullptr);
    }

    Tally tally;
    json rows = json::array();
    json checks = json::array();
    for (std::size_t r = 0; r < kRows.size(); ++r) {
        const RowPlan& p = kRows[r];
        RowRecord rec;
        rec.planned = apps::pilotScaled(p.experiments, o.pilot);
        const std::vector<st::CaseCTruth> truths = rowTruths(p, rec.planned, pool.omega_n, o);
        riptide::random::Engine engine = riptide::random::perEventEngine(o.split_seed, p.section);
        const st::ResampledRow split =
            st::resampleRow(pool.sample.truth_En, g, truths, p.events, engine);
        rec.experiments = split.experiments.size();
        rec.resampling = resamplingJson(split, p.events);
        std::cout << "  " << rowName(p, rec.experiments) << " of " << rec.planned << ", "
                  << split.visited_total << " events visited\n";
        if (setup && rec.experiments >= apps::kMinExperiments) {
            const TwinRow tw = twinRow(twin->at("rows").at(r));
            if (tw.experiments != rec.experiments) {
                throw std::invalid_argument("the twin has " + std::to_string(tw.experiments) +
                                            " experiments in " + rowName(p, rec.experiments));
            }
            const Stopwatch watch(rowName(p, rec.experiments));
            rec.ensemble = st::summariseCaseCEnsemble(
                analyseGeant4Row(*setup, pool, truths, split, o.threads));
            watch.report();
            // N = 1000 is judged on the repeat; the M = 40 row is reported
            // without a verdict that counts.
            const Role role =
                p.events == kBiasEvents && !p.repeat ? Role::reference : Role::counted;
            checks.push_back(geant4Checks(tally, p, *rec.ensemble, tw, pool.omega_n, role));
        } else if (setup) {
            checks.push_back(nullptr);
        }
        rows.push_back(rowJson(p, rec));
    }

    json record = recordHeader(o, caseCDomain());
    record["source"] = pool.source;
    record["source"]["pool_events"] = pool.sample.size();
    record["source"]["pool_density"] = poolDensityJson(g);
    record["truth"] = {{"omega_n", directionJson(pool.omega_n)}};
    record["settings"]["acceptance"] = acceptance.info;
    record["rows"] = rows;
    if (setup) {
        record["twin_id"] = o.twin_id;
        record["checks"] = checks;
        record["tally"] = tally.summary();
        std::cout << "checks: " << tally.counted_failures.size() << " of " << tally.counted
                  << " counted items failed, " << tally.reference_failures.size()
                  << " reference items out\n";
    }
    return record;
}

// The toy twin of a pool record: its truths, N and M per row.
json toyRecord(const fs::path& root, const Options& o) {
    const json match = io::readJsonFile(riptide::data::statsDir(root, o.match_id) / kCaseCName);
    if (match.at("source").at("kind").get<std::string>() != "geant4") {
        throw std::invalid_argument("--match " + o.match_id + " is not a Geant4 record");
    }
    requireSameSettings(match, o, o.match_id);
    const Direction omega_n = directionFromJson(match.at("truth").at("omega_n"));
    const Acceptance acceptance = loadAcceptance(root, o.acceptance_id);
    const st::AcceptanceRatio* ratio = acceptance.ratio ? &*acceptance.ratio : nullptr;
    const Setup setup(o.resolution, ratio);

    json rows = json::array();
    for (std::size_t r = 0; r < kRows.size(); ++r) {
        const RowPlan& p = kRows[r];
        RowRecord rec;
        rec.planned = apps::pilotScaled(p.experiments, o.pilot);
        rec.experiments = match.at("rows").at(r).at("experiments").get<std::size_t>();
        rec.resampling = nullptr;
        if (rec.experiments >= apps::kMinExperiments) {
            const std::vector<st::CaseCTruth> truths = rowTruths(p, rec.experiments, omega_n, o);
            const Stopwatch watch(rowName(p, rec.experiments));
            rec.ensemble = st::summariseCaseCEnsemble(
                analyseToyRow(setup, p, truths, rec.experiments, ratio, o));
            watch.report();
        }
        rows.push_back(rowJson(p, rec));
    }
    json record = recordHeader(o, setup.domain);
    record["source"] = {{"kind", "toy"},
                        {"match",
                         {{"stats_id", o.match_id},
                          {"run_id", match.at("source").at("run_id")},
                          {"mode", match.at("mode")}}}};
    record["truth"] = {{"omega_n", directionJson(omega_n)}};
    record["settings"]["acceptance"] = acceptance.info;
    record["rows"] = rows;
    return record;
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

    std::cout << kAppName << ' ' << riptide::kProjectVersion << " (" << riptide::kGitDescribe
              << "), " << o.threads << " threads, " << modeName(o.mode())
              << (o.pilot ? ", pilot" : "") << '\n';
    const json record = o.mode() == Mode::toy ? toyRecord(root, o) : geant4Record(root, o);
    fs::create_directories(out_dir);
    io::writeFileAtomically(out_dir / kCaseCName, record.dump(2) + "\n");
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
        lyra::opt(o.run_id, "id")["--run-id"]("Geant4 source: the flat riptide_sim run") |
        lyra::opt(o.hp_check_id,
                  "id")["--hp-from-check"]("check with the HP mass (needed by run.json v1)") |
        lyra::opt(o.pool_only)["--pool-only"]("the resampling only, no reconstruction") |
        lyra::opt(o.toy)["--toy"]("the toy twin of the record named by --match") |
        lyra::opt(o.match_id, "id")["--match"]("stats id of the Geant4 record the toy twins") |
        lyra::opt(o.twin_id, "id")["--twin"]("stats id of the toy twin the rms are compared with") |
        lyra::opt(o.seed, "n")["--seed"]("seed of truths, smearing and toy (default 20260907)") |
        lyra::opt(o.split_seed, "n")["--split-seed"]("seed of the resampling (default 20260907)") |
        lyra::opt(o.resolution.sigma_Ep, "MeV")["--sigma-Ep"]("Ep resolution (default 0.10)") |
        lyra::opt(o.resolution.sigma_theta,
                  "rad")["--sigma-theta"]("track resolution per component (default 0.08)") |
        lyra::opt(o.acceptance_id,
                  "id")["--acceptance"]("stats id of acceptance.json: selection weights (D40)") |
        lyra::opt(o.threads, "n")["--threads"]("threads across experiments (default 4)") |
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
        if (o.threads == 0) {
            throw apps::CliError("--threads must be at least 1");
        }
        requireConsistent(o);
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
