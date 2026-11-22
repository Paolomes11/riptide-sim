// Times the Case A likelihood on pseudo-events: one event and a block of
// events, best of a few repetitions, single thread. The numbers belong in the
// validation notes, not in a test: wall times on shared machines are noisy.
// --dump writes the log-likelihood matrix as raw doubles, so two builds can be
// compared number by number after an optimisation.

#include <lyra/lyra.hpp>

#include <algorithm>
#include <chrono>
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
#include "riptide/common/version.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment.hpp"

namespace {

namespace st = riptide::stats;
using Clock = std::chrono::steady_clock;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
constexpr std::uint64_t kProjectSeed = 20260907;
constexpr std::size_t kDefaultEvents = 1000;
constexpr int kDefaultRepeats = 3;
constexpr st::CaseAResolution kToyResolution{0.10, 0.08}; // MeV, rad

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
        std::ofstream out(dump, std::ios::binary);
        out.write(reinterpret_cast<const char*>(matrix.values.data()),
                  static_cast<std::streamsize>(matrix.values.size() * sizeof(double)));
        if (!out) {
            throw std::runtime_error("cannot write " + dump);
        }
    }
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    bool show_help = false;
    std::size_t events = kDefaultEvents;
    int repeats = kDefaultRepeats;
    std::string dump;
    const auto cli =
        lyra::cli() | lyra::help(show_help) |
        lyra::opt(events, "n")["--events"]("events in the timed block (default 1000)") |
        lyra::opt(repeats, "n")["--repeats"]("repetitions, the best is kept (default 3)") |
        lyra::opt(dump, "file")["--dump"]("write the log-likelihood matrix as raw doubles");
    const auto result = cli.parse(lyra::args(argc, argv));
    if (!result || events == 0 || repeats < 1) {
        std::cerr << "riptide_stats_bench: "
                  << (result ? std::string("need --events >= 1 and --repeats >= 1")
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
        return bench(events, repeats, dump);
    } catch (const std::exception& e) {
        std::cerr << "riptide_stats_bench: " << e.what() << '\n';
        return kExitRuntime;
    }
}
