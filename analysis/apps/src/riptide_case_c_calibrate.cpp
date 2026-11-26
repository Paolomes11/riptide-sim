// Calibrates the Case C engine with pseudo-experiments drawn from the toy
// model (D47), so the port can be compared with the toy before any Geant4
// data. It runs:
//   - the checklist on the (N, M) pairs of the toy plus a repeat at N = 1000;
//   - the robustness to the prior uniform in sigma_E at N = 150, on the very
//     experiments of the checklist row;
//   - the edge rows (mu_E near En_max, truth truncated: the model is exact);
//   - the spill rows (truth not truncated: a sensitivity test, expected to
//     fail).
// Every number with a threshold of the written predictions carries its
// verdict next to it, and the verdicts are tallied per group, so the record
// says by itself how many items failed against the false alarms expected.
// Experiments run on several threads; each has its own engines and output
// slot, so the numbers do not depend on --threads. It writes case_c.json
// under the stats directory of the data root; the file holds no wall time.

#include <lyra/lyra.hpp>
#include <nlohmann/json.hpp>

#include <array>
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
#include "riptide/stats/calibration_c.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_c.hpp"
#include "riptide/stats/case_c_posterior.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_c.hpp"

namespace {

namespace fs = std::filesystem;
namespace st = riptide::stats;
namespace units = riptide::units;
using Clock = std::chrono::steady_clock;
using nlohmann::json;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
constexpr std::uint64_t kProjectSeed = 20260907;
constexpr st::CaseBResolution kToyResolution{0.10, 0.08}; // MeV, rad
constexpr unsigned kDefaultThreads = 4;
constexpr int kCalibrationVersion = 1;
constexpr const char* kCalibrationName = "case_c.json";
// Generous upper bound on the output: one small JSON file.
constexpr std::uintmax_t kOutputBytes = std::uintmax_t{1} << 20U; // 1 MiB

// A checklist row of the toy and its reference numbers (the toy run the
// predictions quote): rms of the estimates of mu_E (MeV), log sigma_E and
// omega_n (degrees).
struct ToyRow {
    std::size_t events;      // N per experiment
    std::size_t experiments; // M
    double rms_mu;
    double rms_log_sigma;
    double rms_omega_deg;
};
const std::vector<ToyRow> kToyRows{{50, 200, 0.553, 0.805, 1.85},
                                   {150, 200, 0.350, 0.264, 0.88},
                                   {300, 100, 0.207, 0.113, 0.69},
                                   {1000, 40, 0.092, 0.066, 0.34}};
// The repeat at N = 1000 on disjoint ids: the verdicts at N = 1000 are read
// on it, and its rms is compared with the toy row at M = 40.
constexpr std::size_t kRepeatExperiments = 240;
// The truncation rows of D47, edge and spill alike.
struct SizedRow {
    std::size_t events;
    std::size_t experiments;
};
const std::vector<SizedRow> kTruncationRows{{300, 100}, {1000, 40}};

// Disjoint event-id sections, one per row, so no two draws share an engine.
constexpr std::uint64_t kSectionChecklist = 16;
constexpr std::uint64_t kSectionRepeat = 24;
constexpr std::uint64_t kSectionEdge = 40;
constexpr std::uint64_t kSectionSpill = 48;

// Thresholds of the predictions, section 2.
constexpr double kBiasStdFraction = 0.25;    // |bias| < this x mean posterior std
constexpr double kToleranceSigmas = 2.0;     // statistical bounds at 2 sigma
constexpr double kPullWidthTolerance = 0.07; // for M >= kPullWidthMinExperiments
constexpr std::size_t kPullWidthMinExperiments = 200;
constexpr std::size_t kSlopeMinEvents = 150;
constexpr double kOmegaSlope = -0.5; // 1 / sqrt(N)
constexpr double kOmegaSlopeTolerance = 0.07;
constexpr double kToySlopeMu = -0.71;
constexpr double kToySlopeLogSigma = -0.73;
constexpr double kToySlopeTolerance = 0.15;
// The declared exception: the log sigma_E coverage at N = 50 is low in the
// toy too, a likely prior effect; it is compared with the toy, without a
// verdict that counts.
constexpr std::size_t kLowCoverageEvents = 50;
constexpr std::array<double, 3> kToyLogSigmaCoverageLowN{0.61, 0.81, 0.87};
// Section 3: the prior uniform in sigma_E at N = 150.
constexpr std::size_t kPriorEvents = 150;
constexpr double kPriorMuShiftMax = 0.5;     // posterior stds
constexpr double kPriorLogSigmaShift = 0.21; // posterior stds, toy mean shift
constexpr double kPriorLogSigmaShiftTol = 0.10;
// Section 3: the spill row at N = 1000 breaks the 68% coverage of mu_E as in
// the toy.
constexpr std::size_t kSpillCheckEvents = 1000;
constexpr double kSpillMuCoverageMin = 0.38;
constexpr double kSpillMuCoverageMax = 0.55;
// Section 2: failures expected without a defect, and the count that opens a
// diagnosis.
constexpr std::size_t kExpectedFalseAlarmsMin = 2;
constexpr std::size_t kExpectedFalseAlarmsMax = 4;
constexpr std::size_t kDiagnosisAbove = 7;

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

// The role of a verdict: `counted` items make the false-alarm statistics of
// their group; a `required` item must pass in every run (no statistical
// slack); a `reference` item is a comparison reported without a verdict that
// counts.
enum class Role { counted, required, reference };

const char* roleName(Role r) {
    switch (r) {
    case Role::counted:
        return "counted";
    case Role::required:
        return "required";
    case Role::reference:
        return "reference";
    }
    return "reference";
}

struct Tally {
    std::size_t counted = 0;
    std::vector<std::string> counted_failures;
    std::vector<std::string> required_failures;
    std::vector<std::string> reference_failures;

    json check(const std::string& name, double value, double lo, double hi, Role role) {
        const bool pass = value >= lo && value <= hi;
        if (role == Role::counted) {
            ++counted;
        }
        if (!pass) {
            switch (role) {
            case Role::counted:
                counted_failures.push_back(name);
                break;
            case Role::required:
                required_failures.push_back(name);
                break;
            case Role::reference:
                reference_failures.push_back(name);
                break;
            }
        }
        return {{"value", value}, {"lo", lo}, {"hi", hi}, {"pass", pass}, {"role", roleName(role)}};
    }

    [[nodiscard]] json summary() const {
        return {{"counted", counted},
                {"counted_failures", counted_failures},
                {"required_failures", required_failures},
                {"reference_failures", reference_failures}};
    }
};

double binomialBound(double p, std::size_t n) {
    return kToleranceSigmas * std::sqrt(p * (1.0 - p) / static_cast<double>(n));
}

// Relative error of an rms over M experiments is about 1 / sqrt(2M); the
// ratio of two independent rms values adds the two in quadrature.
double rmsRatioBound(std::size_t m, std::size_t m_toy) {
    return kToleranceSigmas *
           std::sqrt((0.5 / static_cast<double>(m)) + (0.5 / static_cast<double>(m_toy)));
}

// The width of the pulls has a relative error of about 1 / sqrt(2M); from
// M = 200 the predictions fix the bound at 0.07.
double pullWidthBound(std::size_t m) {
    if (m >= kPullWidthMinExperiments) {
        return kPullWidthTolerance;
    }
    return kToleranceSigmas / std::sqrt(2.0 * static_cast<double>(m));
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

json directionJson(const st::CaseBEnsemble& e, const st::DirectionBias& b) {
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
            {"tangent_bias_error_rad", b.error}};
}

json ensembleJson(const st::CaseCEnsemble& e) {
    return {{"mu_MeV", hyperJson(e.mu)},
            {"log_sigma", hyperJson(e.log_sigma)},
            {"omega", directionJson(e.direction, e.direction_bias)}};
}

struct Setup {
    st::EnergyDomain domain;
    st::CaseCModel model;
    std::vector<riptide::geometry::Direction> coarse_grid = st::sphereGrid(st::kSphereGridPixels);

    explicit Setup(const st::CaseBResolution& r)
        : model(st::energyGrid(domain), st::recoilAngleGrid(st::kTrackAngleGridPoints),
                st::Grid{0.0, units::pi, st::kTrackObservedGridPoints}, r) {}
};

struct RowResult {
    st::CaseCRowSpec spec;
    std::vector<st::CaseCExperimentResult> results;
    st::CaseCEnsemble ensemble;
};

RowResult runRow(const Setup& s, const Options& o, st::CaseCRowSpec spec, const std::string& name) {
    spec.experiments = riptide::apps::pilotScaled(spec.experiments, o.pilot);
    const Stopwatch watch(name + " N = " + std::to_string(spec.events) +
                          ", M = " + std::to_string(spec.experiments));
    RowResult out{
        spec,
        st::runCaseCRow(s.model, s.domain, spec, o.resolution, o.seed, s.coarse_grid, o.threads),
        {}};
    out.ensemble = st::summariseCaseCEnsemble(out.results);
    watch.report(spec.experiments);
    return out;
}

std::string label(const RowResult& r, const std::string& what) {
    return "N=" + std::to_string(r.spec.events) + " M=" + std::to_string(r.spec.experiments) + " " +
           what;
}

const std::array<const char*, 3> kLevelNames{"0.68", "0.90", "0.95"};

// Pull widths and coverages of the three parameters, the items the edge rows
// share with the checklist.
json calibrationChecks(Tally& t, const RowResult& r, Role role) {
    const st::CaseCEnsemble& e = r.ensemble;
    const std::size_t m = r.spec.experiments;
    const double pull_bound = pullWidthBound(m);
    json out;
    out["pull_width_mu"] = t.check(label(r, "pull width mu"), e.mu.pull_width, 1.0 - pull_bound,
                                   1.0 + pull_bound, role);
    out["pull_width_log_sigma"] = t.check(label(r, "pull width log sigma"), e.log_sigma.pull_width,
                                          1.0 - pull_bound, 1.0 + pull_bound, role);
    out["pull_rms_omega"] = t.check(label(r, "pull rms omega"), e.direction.pull_rms,
                                    1.0 - pull_bound, 1.0 + pull_bound, role);
    const bool low_n = r.spec.events == kLowCoverageEvents;
    json coverage = json::array();
    for (std::size_t l = 0; l < st::kHyperCredibleLevels.size(); ++l) {
        const double p = st::kHyperCredibleLevels[l];
        const double bound = binomialBound(p, m);
        const std::string level = kLevelNames[l];
        json c;
        c["level"] = p;
        c["mu"] =
            t.check(label(r, "coverage mu " + level), e.mu.coverage[l], p - bound, p + bound, role);
        // The declared exception: against the toy, not the nominal level.
        const double ls_ref = low_n ? kToyLogSigmaCoverageLowN[l] : p;
        const double ls_bound = binomialBound(ls_ref, m);
        c["log_sigma"] =
            t.check(label(r, "coverage log sigma " + level), e.log_sigma.coverage[l],
                    ls_ref - ls_bound, ls_ref + ls_bound, low_n ? Role::reference : role);
        c["omega"] = t.check(label(r, "coverage omega " + level), e.direction.coverage[l],
                             p - bound, p + bound, role);
        coverage.push_back(c);
    }
    out["coverage"] = coverage;
    out["truth_in_cap"] = t.check(label(r, "truth in cap"), e.direction.truth_in_cap_fraction, 1.0,
                                  1.0, role == Role::reference ? role : Role::required);
    return out;
}

// All the section-2 items of one row against its toy row.
json checklistChecks(Tally& t, const RowResult& r, const ToyRow& toy, Role role) {
    const st::CaseCEnsemble& e = r.ensemble;
    json out = calibrationChecks(t, r, role);
    const double mu_bound = kBiasStdFraction * e.mu.std_mean;
    const double ls_bound = kBiasStdFraction * e.log_sigma.std_mean;
    out["bias_mu"] = t.check(label(r, "bias mu"), e.mu.bias, -mu_bound, mu_bound, role);
    out["bias_log_sigma"] =
        t.check(label(r, "bias log sigma"), e.log_sigma.bias, -ls_bound, ls_bound, role);
    const double ratio_bound = rmsRatioBound(r.spec.experiments, toy.experiments);
    const double lo = 1.0 - ratio_bound;
    const double hi = 1.0 + ratio_bound;
    out["toy_events"] = toy.events;
    out["toy_experiments"] = toy.experiments;
    out["rms_ratio_mu"] = t.check(label(r, "rms ratio mu"), e.mu.rms / toy.rms_mu, lo, hi, role);
    out["rms_ratio_log_sigma"] =
        t.check(label(r, "rms ratio log sigma"), e.log_sigma.rms / toy.rms_log_sigma, lo, hi, role);
    out["rms_ratio_omega"] =
        t.check(label(r, "rms ratio omega"),
                e.direction.error_rms / (toy.rms_omega_deg * units::deg), lo, hi, role);
    return out;
}

json rowJson(const RowResult& r, json checks) {
    const st::CaseCRowSpec& s = r.spec;
    return {{"events", s.events},
            {"experiments", s.experiments},
            {"section", s.section},
            {"truth",
             {{"mu_min_MeV", s.truth.mu_min},
              {"mu_max_MeV", s.truth.mu_max},
              {"sigma_min_MeV", s.truth.sigma_min},
              {"sigma_max_MeV", s.truth.sigma_max}}},
            {"spectrum", s.spectrum == st::CaseCSpectrum::truncated ? "truncated" : "untruncated"},
            {"ensemble", ensembleJson(r.ensemble)},
            {"checks", std::move(checks)}};
}

struct Slopes {
    double mu;
    double log_sigma;
    double omega;
};

Slopes slopesOf(const std::vector<const RowResult*>& rows) {
    std::vector<double> events;
    std::vector<double> mu;
    std::vector<double> log_sigma;
    std::vector<double> omega;
    for (const RowResult* r : rows) {
        if (r->spec.events >= kSlopeMinEvents) {
            events.push_back(static_cast<double>(r->spec.events));
            mu.push_back(r->ensemble.mu.rms);
            log_sigma.push_back(r->ensemble.log_sigma.rms);
            omega.push_back(r->ensemble.direction.error_rms);
        }
    }
    return {st::logLogSlope(events, mu), st::logLogSlope(events, log_sigma),
            st::logLogSlope(events, omega)};
}

json slopeChecks(Tally& t, const Slopes& s, const std::string& name, Role role) {
    return {{"mu", t.check(name + " slope mu", s.mu, kToySlopeMu - kToySlopeTolerance,
                           kToySlopeMu + kToySlopeTolerance, role)},
            {"log_sigma",
             t.check(name + " slope log sigma", s.log_sigma, kToySlopeLogSigma - kToySlopeTolerance,
                     kToySlopeLogSigma + kToySlopeTolerance, role)},
            {"omega", t.check(name + " slope omega", s.omega, kOmegaSlope - kOmegaSlopeTolerance,
                              kOmegaSlope + kOmegaSlopeTolerance, role)}};
}

// The common mean pull of mu_E over the rows, weighted by their experiments;
// the predictions report it without a threshold.
json commonMuPull(const std::vector<const RowResult*>& rows) {
    double weighted = 0.0;
    double variance = 0.0;
    double total = 0.0;
    for (const RowResult* r : rows) {
        const auto m = static_cast<double>(r->spec.experiments);
        weighted += m * r->ensemble.mu.pull_mean;
        variance += m * m * r->ensemble.mu.pull_mean_error * r->ensemble.mu.pull_mean_error;
        total += m;
    }
    return {{"mean", weighted / total}, {"error", std::sqrt(variance) / total}};
}

json priorJson(Tally& t, const RowResult& ref) {
    const st::CaseCEnsemble alt = st::summariseCaseCEnsemble(ref.results, true);
    const st::CaseCPriorShift shift = st::caseCPriorShift(ref.results);
    const std::size_t m = ref.spec.experiments;
    json coverage = json::array();
    for (std::size_t l = 0; l < st::kHyperCredibleLevels.size(); ++l) {
        const double p = st::kHyperCredibleLevels[l];
        const double bound = binomialBound(p, m);
        coverage.push_back(
            t.check(label(ref, "uniform sigma coverage log sigma " + std::string(kLevelNames[l])),
                    alt.log_sigma.coverage[l], p - bound, p + bound, Role::counted));
    }
    const auto shiftJson = [](const st::HyperShift& s) {
        return json{{"mean", s.mean}, {"mean_error", s.mean_error}, {"max_abs", s.max_abs}};
    };
    return {{"events", ref.spec.events},
            {"experiments", m},
            {"prior", "uniform in sigma_E"},
            {"baseline", "uniform in log sigma_E"},
            {"shift_mu", shiftJson(shift.mu)},
            {"shift_log_sigma", shiftJson(shift.log_sigma)},
            {"ensemble", ensembleJson(alt)},
            {"checks",
             {{"shift_mu_max", t.check(label(ref, "uniform sigma max shift mu"), shift.mu.max_abs,
                                       0.0, kPriorMuShiftMax, Role::counted)},
              {"shift_log_sigma_mean",
               t.check(label(ref, "uniform sigma mean shift log sigma"), shift.log_sigma.mean,
                       kPriorLogSigmaShift - kPriorLogSigmaShiftTol,
                       kPriorLogSigmaShift + kPriorLogSigmaShiftTol, Role::counted)},
              {"coverage_log_sigma", coverage}}}};
}

std::vector<RowResult> runTruncationRows(const Setup& s, const Options& o,
                                         const st::CaseCTruthRange& truth,
                                         st::CaseCSpectrum spectrum, std::uint64_t section,
                                         const std::string& name) {
    std::vector<RowResult> rows;
    for (std::size_t r = 0; r < kTruncationRows.size(); ++r) {
        const st::CaseCRowSpec spec{kTruncationRows[r].events,
                                    kTruncationRows[r].experiments,
                                    section + r,
                                    truth,
                                    spectrum,
                                    false};
        rows.push_back(runRow(s, o, spec, name));
    }
    return rows;
}

json tallyJson(const Tally& t) {
    json out = t.summary();
    out["expected_false_alarms"] = {kExpectedFalseAlarmsMin, kExpectedFalseAlarmsMax};
    out["diagnosis_above"] = kDiagnosisAbove;
    return out;
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

    std::cout << "riptide_case_c_calibrate " << riptide::kProjectVersion << " ("
              << riptide::kGitDescribe << "), " << o.threads << " threads"
              << (o.pilot ? ", pilot" : "") << '\n';
    const Setup setup(o.resolution);

    // The checklist; the row at kPriorEvents also carries the alternative
    // prior, on the same base.
    std::vector<RowResult> checklist;
    for (std::size_t r = 0; r < kToyRows.size(); ++r) {
        const ToyRow& toy = kToyRows[r];
        const st::CaseCRowSpec spec{toy.events,
                                    toy.experiments,
                                    kSectionChecklist + r,
                                    st::kCaseCChecklistTruth,
                                    st::CaseCSpectrum::truncated,
                                    toy.events == kPriorEvents};
        checklist.push_back(runRow(setup, o, spec, "checklist"));
    }
    const ToyRow& toy_last = kToyRows.back();
    const RowResult repeat = runRow(setup, o,
                                    {toy_last.events, kRepeatExperiments, kSectionRepeat,
                                     st::kCaseCChecklistTruth, st::CaseCSpectrum::truncated, false},
                                    "repeat");
    const std::vector<RowResult> edge = runTruncationRows(
        setup, o, st::kCaseCEdgeTruth, st::CaseCSpectrum::truncated, kSectionEdge, "edge");
    const std::vector<RowResult> spill = runTruncationRows(
        setup, o, st::kCaseCChecklistTruth, st::CaseCSpectrum::untruncated, kSectionSpill, "spill");

    // Section 2: the verdicts at N = 1000 are read on the repeat; the M = 40
    // row is the reference of the comparison with the toy.
    Tally port;
    json rows_json = json::array();
    std::vector<const RowResult*> rows;
    std::vector<const RowResult*> rows_with_repeat;
    for (std::size_t r = 0; r < checklist.size(); ++r) {
        const bool last = r + 1 == checklist.size();
        rows_json.push_back(
            rowJson(checklist[r], checklistChecks(port, checklist[r], kToyRows[r],
                                                  last ? Role::reference : Role::counted)));
        rows.push_back(&checklist[r]);
        rows_with_repeat.push_back(last ? &repeat : &checklist[r]);
    }
    const json repeat_json =
        rowJson(repeat, checklistChecks(port, repeat, toy_last, Role::counted));
    json contraction = {{"slope_min_events", kSlopeMinEvents},
                        {"slopes_with_repeat", slopeChecks(port, slopesOf(rows_with_repeat),
                                                           "with repeat", Role::counted)},
                        {"slopes", slopeChecks(port, slopesOf(rows), "toy rows", Role::reference)}};

    Tally prior;
    const RowResult* prior_row = nullptr;
    for (const RowResult& r : checklist) {
        if (r.spec.with_alternate) {
            prior_row = &r;
        }
    }
    const json prior_json = priorJson(prior, *prior_row);

    Tally edge_tally;
    json edge_json = json::array();
    for (const RowResult& r : edge) {
        edge_json.push_back(rowJson(r, calibrationChecks(edge_tally, r, Role::counted)));
    }
    // The spill rows are a sensitivity test: everything is reference.
    Tally spill_tally;
    json spill_json = json::array();
    for (const RowResult& r : spill) {
        json checks = calibrationChecks(spill_tally, r, Role::reference);
        if (r.spec.events == kSpillCheckEvents) {
            checks["toy_coverage_mu_0.68"] =
                spill_tally.check(label(r, "toy range coverage mu 0.68"), r.ensemble.mu.coverage[0],
                                  kSpillMuCoverageMin, kSpillMuCoverageMax, Role::reference);
        }
        spill_json.push_back(rowJson(r, std::move(checks)));
    }

    json record;
    record["calibration_version"] = kCalibrationVersion;
    record["calibration_id"] = o.calibration_id;
    record["project_version"] = riptide::kProjectVersion;
    record["git_describe"] = riptide::kGitDescribe;
    // The thread count is left out: it does not change a number.
    record["settings"] = {{"seed", o.seed},
                          {"sigma_Ep_MeV", o.resolution.sigma_Ep},
                          {"sigma_theta_rad", o.resolution.sigma_theta},
                          {"pilot", o.pilot},
                          {"levels", st::kHyperCredibleLevels},
                          {"En_domain", domainJson(setup.domain)},
                          {"theta_points", st::kTrackAngleGridPoints},
                          {"theta_obs_points", st::kTrackObservedGridPoints},
                          {"sphere_pixels", st::kSphereGridPixels},
                          {"cap_pixels", st::kCapGridPixels},
                          {"cap_delta_log", st::kCapDeltaLog},
                          {"hyper_global_points", {st::kHyperMuPoints, st::kHyperSigmaPoints}},
                          {"hyper_sigma_range_MeV", {st::kHyperSigmaMin, st::kHyperSigmaMax}},
                          {"hyper_fine_points", st::kHyperFinePoints},
                          {"hyper_window_delta_log", st::kHyperWindowDeltaLog}};
    record["checklist"] = {{"rows", rows_json},
                           {"repeat", repeat_json},
                           {"contraction", contraction},
                           {"common_mu_pull", commonMuPull(rows_with_repeat)},
                           {"tally", tallyJson(port)}};
    record["prior"] = {{"row", prior_json}, {"tally", prior.summary()}};
    record["edge"] = {{"rows", edge_json}, {"tally", edge_tally.summary()}};
    record["spill"] = {{"rows", spill_json}, {"tally", spill_tally.summary()}};

    fs::create_directories(out_dir);
    riptide::io::writeFileAtomically(out_dir / kCalibrationName, record.dump(2) + "\n");
    std::cout << "checklist: " << port.counted_failures.size() << " of " << port.counted
              << " counted items failed, " << port.required_failures.size() << " required failed\n"
              << "prior: " << prior.counted_failures.size() << " of " << prior.counted
              << "; edge: " << edge_tally.counted_failures.size() << " of " << edge_tally.counted
              << ", " << edge_tally.required_failures.size() << " required failed\n"
              << "riptide_case_c_calibrate: output in " << out_dir << '\n';
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
        std::cerr << "riptide_case_c_calibrate: " << result.message() << "\n\n" << cli << '\n';
        return kExitUsage;
    }
    if (o.threads == 0) {
        std::cerr << "riptide_case_c_calibrate: --threads must be at least 1\n";
        return kExitUsage;
    }
    try {
        return calibrate(o);
    } catch (const std::exception& e) {
        std::cerr << "riptide_case_c_calibrate: " << e.what() << '\n';
        return kExitRuntime;
    }
}
