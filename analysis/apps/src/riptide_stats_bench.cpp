// Times the Case A likelihood on pseudo-events (one event and a block of
// events), one Case B experiment phase by phase (track-energy table, coarse
// sphere plus cap, summaries), or one Case C experiment stage by stage
// (stage 1, per-event base, stage 2, summaries and predictive, stage 1');
// best of a few repetitions, single thread. The numbers belong in the validation notes, not in a
// test: wall times on shared machines are noisy. --dump writes the computed numbers as raw doubles,
// so two builds can be compared number by number after an optimisation: the log-likelihood matrix
// for Case A; the table, the cap log-posterior and omega_hat for Case B; the base, the fine
// hyperparameter log-posterior, the predictive moments and omega_hat for Case C.

#include <lyra/lyra.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "riptide/common/git_info.hpp"
#include "riptide/common/units.hpp"
#include "riptide/common/version.hpp"
#include "riptide/stats/calibration.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_b_posterior.hpp"
#include "riptide/stats/case_c.hpp"
#include "riptide/stats/case_c_posterior.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"
#include "riptide/stats/pseudo_experiment_c.hpp"

namespace {

namespace st = riptide::stats;
using Clock = std::chrono::steady_clock;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
constexpr std::uint64_t kProjectSeed = 20260907;
constexpr std::size_t kDefaultEvents = 1000;
constexpr int kDefaultRepeats = 3;
constexpr st::CaseAResolution kToyResolution{0.10, 0.08};  // MeV, rad
constexpr st::CaseBResolution kToyResolutionB{0.10, 0.08}; // MeV, rad

struct Timing {
    double best_s = std::numeric_limits<double>::infinity();
    double checksum = 0.0; // keeps the work observable and flags a changed result
};

Timing timeBlock(const st::CaseAModel& model, const st::PseudoSample& sample, std::size_t n,
                 int repeats, st::LogLikelihoodMatrix* keep) {
    const auto count = static_cast<std::ptrdiff_t>(n);
    const std::vector<double> Ep(sample.reco_Ep.begin(), sample.reco_Ep.begin() + count);
    const std::vector<double> theta(sample.reco_theta_p.begin(),
                                    sample.reco_theta_p.begin() + count);
    Timing timing;
    for (int r = 0; r < repeats; ++r) {
        const auto start = Clock::now();
        st::LogLikelihoodMatrix matrix = model.logLikelihoods(Ep, theta);
        const std::chrono::duration<double> elapsed = Clock::now() - start;
        timing.best_s = std::min(timing.best_s, elapsed.count());
        timing.checksum = 0.0;
        for (const double v : matrix.values) {
            timing.checksum += v;
        }
        if (keep != nullptr) {
            *keep = std::move(matrix);
        }
    }
    return timing;
}

void writeDump(const std::string& dump, const std::vector<double>& values) {
    std::ofstream out(dump, std::ios::binary);
    out.write(reinterpret_cast<const char*>(values.data()),
              static_cast<std::streamsize>(values.size() * sizeof(double)));
    if (!out) {
        throw std::runtime_error("cannot write " + dump);
    }
}

double elapsedSince(Clock::time_point start) {
    const std::chrono::duration<double> elapsed = Clock::now() - start;
    return elapsed.count();
}

// Best time of each phase over the repetitions, taken separately.
struct CaseBTiming {
    double table_s = std::numeric_limits<double>::infinity();
    double posterior_s = std::numeric_limits<double>::infinity();
    double summary_s = std::numeric_limits<double>::infinity();
    double total_s = std::numeric_limits<double>::infinity();
};

int benchCaseB(std::size_t events, int repeats, const std::string& dump) {
    const st::EnergyDomain domain;
    const st::Grid En = st::energyGrid(domain);
    const st::CaseBModel model(En, st::recoilAngleGrid(st::kTrackAngleGridPoints),
                               st::Grid{0.0, riptide::units::pi, st::kTrackObservedGridPoints},
                               st::flatEnergyLogPrior(En), kToyResolutionB);
    const std::vector<riptide::geometry::Direction> coarse = st::sphereGrid(st::kSphereGridPixels);
    const st::CaseBPseudoExperiment experiment =
        st::caseBPseudoExperiment(domain, events, kToyResolutionB, kProjectSeed, 0, 0);
    const st::DirectionLogPrior uniform = st::uniformDirectionLogPrior;

    CaseBTiming timing;
    st::TrackEnergyTable table;
    st::DirectionPosterior posterior;
    st::DirectionSummary summary;
    for (int r = 0; r < repeats; ++r) {
        const auto start = Clock::now();
        table = model.trackEnergyTable(experiment.sample.reco_Ep);
        const double table_s = elapsedSince(start);
        const auto posterior_start = Clock::now();
        posterior =
            st::refineSharedDirection(model, table, experiment.sample.reco_track, coarse, uniform);
        const double posterior_s = elapsedSince(posterior_start);
        const auto summary_start = Clock::now();
        summary = st::summarizeDirection(posterior, experiment.truth_omega_n);
        const double summary_s = elapsedSince(summary_start);
        const double total_s = elapsedSince(start);
        timing.table_s = std::min(timing.table_s, table_s);
        timing.posterior_s = std::min(timing.posterior_s, posterior_s);
        timing.summary_s = std::min(timing.summary_s, summary_s);
        timing.total_s = std::min(timing.total_s, total_s);
    }
    double table_sum = 0.0;
    for (const double v : table.values) {
        table_sum += v;
    }
    // -inf pixels (no compatible direction) would swamp the sum.
    double cap_sum = 0.0;
    for (const double v : posterior.cap_log_posterior) {
        cap_sum += std::isfinite(v) ? v : 0.0;
    }

    std::cout.precision(17);
    std::cout << "riptide_stats_bench " << riptide::kProjectVersion << " (" << riptide::kGitDescribe
              << ")\n"
              << "case b, " << events << " events, grid " << En.size() << " x "
              << model.recoilAngleGrid().size() << " x " << model.observedAngleGrid().size()
              << ", sphere " << coarse.size() << ", cap " << posterior.cap.size() << ", best of "
              << repeats << '\n'
              << "table:      " << timing.table_s << " s\n"
              << "posterior:  " << timing.posterior_s << " s\n"
              << "summary:    " << timing.summary_s << " s\n"
              << "total:      " << timing.total_s << " s\n"
              << "checksum table " << table_sum << ", cap " << cap_sum << '\n'
              << "omega_hat " << posterior.omega_hat.x << ' ' << posterior.omega_hat.y << ' '
              << posterior.omega_hat.z << ", error " << summary.error << " rad, sigma "
              << summary.sigma << " rad\n";
    if (!dump.empty()) {
        std::vector<double> values = table.values;
        values.insert(values.end(), posterior.cap_log_posterior.begin(),
                      posterior.cap_log_posterior.end());
        values.push_back(posterior.omega_hat.x);
        values.push_back(posterior.omega_hat.y);
        values.push_back(posterior.omega_hat.z);
        writeDump(dump, values);
    }
    return kExitOk;
}

// Best time of each stage over the repetitions, taken separately.
struct CaseCTiming {
    double stage1_s = std::numeric_limits<double>::infinity();
    double base_s = std::numeric_limits<double>::infinity();
    double stage2_s = std::numeric_limits<double>::infinity();
    double summary_s = std::numeric_limits<double>::infinity();
    double stage1p_s = std::numeric_limits<double>::infinity();
    double total_s = std::numeric_limits<double>::infinity();
};

double finiteSum(const std::vector<double>& values) {
    // -inf entries (no compatible configuration) would swamp the sum.
    double sum = 0.0;
    for (const double v : values) {
        sum += std::isfinite(v) ? v : 0.0;
    }
    return sum;
}

int benchCaseC(std::size_t events, int repeats, const std::string& dump) {
    const st::EnergyDomain domain;
    const st::CaseCModel model(
        st::energyGrid(domain), st::recoilAngleGrid(st::kTrackAngleGridPoints),
        st::Grid{0.0, riptide::units::pi, st::kTrackObservedGridPoints}, kToyResolutionB);
    const std::vector<riptide::geometry::Direction> coarse = st::sphereGrid(st::kSphereGridPixels);
    const st::CaseCPseudoExperiment experiment =
        st::caseCPseudoExperiment(domain, events, kToyResolutionB, st::kCaseCChecklistTruth,
                                  st::CaseCSpectrum::truncated, kProjectSeed, 0, 0);
    const std::vector<double>& Ep = experiment.sample.reco_Ep;
    const std::vector<riptide::geometry::Direction>& tracks = experiment.sample.reco_track;
    const st::DirectionLogPrior uniform = st::uniformDirectionLogPrior;
    const st::HyperLogPrior baseline = st::uniformLogSigmaHyperPrior;

    CaseCTiming timing;
    st::DirectionPosterior stage1;
    st::HierarchicalBase base;
    // Not default-constructible: its grids are built by the model.
    std::optional<st::HyperPosterior> hyper;
    st::CaseCHyperSummary summary;
    st::PredictiveEnergy predictive;
    st::DirectionPosterior stage1p;
    for (int r = 0; r < repeats; ++r) {
        const auto start = Clock::now();
        const st::TrackEnergyTable table = model.stage1Model().trackEnergyTable(Ep);
        stage1 = st::refineSharedDirection(model.stage1Model(), table, tracks, coarse, uniform);
        const double stage1_s = elapsedSince(start);
        const auto base_start = Clock::now();
        base = model.base(Ep, tracks, stage1.omega_hat);
        const double base_s = elapsedSince(base_start);
        const auto stage2_start = Clock::now();
        hyper = model.refineHyperparameters(base, baseline);
        const double stage2_s = elapsedSince(stage2_start);
        const auto summary_start = Clock::now();
        summary = st::summarizeHyperparameters(*hyper);
        predictive = st::predictiveEnergyMoments(*hyper);
        const double summary_s = elapsedSince(summary_start);
        const auto stage1p_start = Clock::now();
        stage1p = st::refineHierarchicalDirection(model, Ep, tracks, coarse, uniform, predictive);
        const double stage1p_s = elapsedSince(stage1p_start);
        const double total_s = elapsedSince(start);
        timing.stage1_s = std::min(timing.stage1_s, stage1_s);
        timing.base_s = std::min(timing.base_s, base_s);
        timing.stage2_s = std::min(timing.stage2_s, stage2_s);
        timing.summary_s = std::min(timing.summary_s, summary_s);
        timing.stage1p_s = std::min(timing.stage1p_s, stage1p_s);
        timing.total_s = std::min(timing.total_s, total_s);
    }
    const st::DirectionSummary direction =
        st::summarizeDirection(stage1p, experiment.truth.omega_n);

    std::cout.precision(17);
    std::cout << "riptide_stats_bench " << riptide::kProjectVersion << " (" << riptide::kGitDescribe
              << ")\n"
              << "case c, " << events << " events, truth mu " << experiment.truth.mu
              << " MeV, sigma " << experiment.truth.sigma << " MeV, hyper grid "
              << st::kHyperMuPoints << " x " << st::kHyperSigmaPoints << " + "
              << st::kHyperFinePoints << " x " << st::kHyperFinePoints << ", best of " << repeats
              << '\n'
              << "stage 1:    " << timing.stage1_s << " s\n"
              << "base:       " << timing.base_s << " s\n"
              << "stage 2:    " << timing.stage2_s << " s\n"
              << "summary:    " << timing.summary_s << " s\n"
              << "stage 1':   " << timing.stage1p_s << " s\n"
              << "total:      " << timing.total_s << " s\n"
              << "checksum base " << finiteSum(base.values) << ", hyper "
              << finiteSum(hyper->log_posterior) << '\n'
              << "mu " << summary.mu.mean << " +- " << summary.mu.std << " MeV, log sigma median "
              << summary.log_sigma.median << ", predictive " << predictive.mean << " +- "
              << predictive.std << " MeV\n"
              << "omega_hat " << stage1p.omega_hat.x << ' ' << stage1p.omega_hat.y << ' '
              << stage1p.omega_hat.z << ", error " << direction.error << " rad, sigma "
              << direction.sigma << " rad\n";
    if (!dump.empty()) {
        std::vector<double> values = base.values;
        values.insert(values.end(), hyper->log_posterior.begin(), hyper->log_posterior.end());
        values.push_back(predictive.mean);
        values.push_back(predictive.std);
        values.push_back(stage1p.omega_hat.x);
        values.push_back(stage1p.omega_hat.y);
        values.push_back(stage1p.omega_hat.z);
        writeDump(dump, values);
    }
    return kExitOk;
}

int bench(std::size_t events, int repeats, const std::string& dump) {
    const st::EnergyDomain domain;
    const st::Grid En = st::energyGrid(domain);
    const st::Grid theta = st::recoilAngleGrid();
    const st::CaseAModel model(En, theta, st::physicalRecoilAngleLogPrior(theta), kToyResolution);
    const st::PseudoSample sample =
        st::pseudoSampleFromPrior(domain, events, kToyResolution, kProjectSeed);

    const Timing one = timeBlock(model, sample, 1, repeats, nullptr);
    st::LogLikelihoodMatrix matrix;
    const Timing block = timeBlock(model, sample, events, repeats, &matrix);

    std::cout.precision(17);
    std::cout << "riptide_stats_bench " << riptide::kProjectVersion << " (" << riptide::kGitDescribe
              << ")\n"
              << "grid " << En.size() << " x " << theta.size() << ", best of " << repeats << '\n'
              << "1 event:      " << one.best_s << " s\n"
              << events << " events: " << block.best_s << " s, checksum " << block.checksum << '\n';
    if (!dump.empty()) {
        writeDump(dump, matrix.values);
    }
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    bool show_help = false;
    std::size_t events = kDefaultEvents;
    int repeats = kDefaultRepeats;
    std::string dump;
    std::string bench_case = "a";
    const auto cli =
        lyra::cli() | lyra::help(show_help) |
        lyra::opt(events, "n")["--events"]("events in the timed block (default 1000)") |
        lyra::opt(repeats, "n")["--repeats"]("repetitions, the best is kept (default 3)") |
        lyra::opt(dump, "file")["--dump"]("write the computed numbers as raw doubles") |
        lyra::opt(bench_case,
                  "a|b|c")["--case"]("Case A likelihood, or one Case B or Case C experiment")
            .choices("a", "b", "c");
    const auto result = cli.parse(lyra::args(argc, argv));
    const bool case_a = bench_case == "a";
    const bool events_ok = events >= 1 && (case_a || events <= st::kMaxExperimentEvents);
    if (!result || !events_ok || repeats < 1) {
        std::cerr << "riptide_stats_bench: "
                  << (result ? std::string(
                                   "need --events in [1, kMaxExperimentEvents] for cases b and c, "
                                   ">= 1 for case a, and --repeats >= 1")
                             : result.message())
                  << "\n\n"
                  << cli << '\n';
        return kExitUsage;
    }
    if (show_help) {
        std::cout << cli << '\n';
        return kExitOk;
    }
    try {
        if (bench_case == "c") {
            return benchCaseC(events, repeats, dump);
        }
        return bench_case == "b" ? benchCaseB(events, repeats, dump) : bench(events, repeats, dump);
    } catch (const std::exception& e) {
        std::cerr << "riptide_stats_bench: " << e.what() << '\n';
        return kExitRuntime;
    }
}
