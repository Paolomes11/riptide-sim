// Times the Case A likelihood on pseudo-events (one event and a block of
// events), or one Case B experiment phase by phase (track-energy table,
// coarse sphere plus cap, summaries); best of a few repetitions, single
// thread. The numbers belong in the validation notes, not in a test: wall
// times on shared machines are noisy. --dump writes the computed numbers as
// raw doubles, so two builds can be compared number by number after an
// optimisation: the log-likelihood matrix for Case A; the table, the cap
// log-posterior and omega_hat for Case B.

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
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"

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
        lyra::opt(bench_case, "a|b")["--case"]("Case A likelihood or one Case B experiment")
            .choices("a", "b");
    const auto result = cli.parse(lyra::args(argc, argv));
    const bool case_b = bench_case == "b";
    const bool events_ok = events >= 1 && (!case_b || events <= st::kMaxExperimentEvents);
    if (!result || !events_ok || repeats < 1) {
        std::cerr << "riptide_stats_bench: "
                  << (result ? std::string("need --events in [1, kMaxExperimentEvents] for case b, "
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
        return case_b ? benchCaseB(events, repeats, dump) : bench(events, repeats, dump);
    } catch (const std::exception& e) {
        std::cerr << "riptide_stats_bench: " << e.what() << '\n';
        return kExitRuntime;
    }
}
