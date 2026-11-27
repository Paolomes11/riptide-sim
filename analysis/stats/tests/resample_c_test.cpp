#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/common/random.hpp"
#include "riptide/stats/acceptance.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_c.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"
#include "riptide/stats/pseudo_experiment_c.hpp"
#include "riptide/stats/resample_c.hpp"

namespace st = riptide::stats;
namespace geo = riptide::geometry;
namespace rnd = riptide::random;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace {

constexpr std::uint64_t kSeed = 20260907;
constexpr double kPi = 3.14159265358979323846;
// The Geant4 domain of D48.
const st::EnergyDomain kDomain{10.0, 19.5, st::kEnergyGridPoints};
const st::CaseBResolution kResolution{0.10, 0.08};
const geo::Direction kOmega{0.0, 0.0, 1.0};

// A pool that is not flat: density rising linearly across the domain, as
// a cross section times an acceptance might shape it, but never zero, so
// every bin is filled. With x = (E - lo) / w the density is (1 + 3x) / 2.5,
// whose inverse CDF is x = (sqrt(1 + 15u) - 1) / 3.
std::vector<double> risingPool(std::size_t n, std::uint64_t seed) {
    auto engine = rnd::perEventEngine(seed, 0);
    const double width = kDomain.En_max - kDomain.En_min;
    std::vector<double> pool(n);
    for (double& En : pool) {
        const double u = rnd::uniform53(engine);
        En = kDomain.En_min + (width * (std::sqrt(1.0 + (15.0 * u)) - 1.0) / 3.0);
    }
    return pool;
}

std::vector<st::CaseCTruth> sameTruths(std::size_t m, double mu, double sigma) {
    st::CaseCTruth truth;
    truth.mu = mu;
    truth.sigma = sigma;
    truth.omega_n = kOmega;
    return std::vector<st::CaseCTruth>(m, truth);
}

// chi^2 of the kept energies on 1 MeV-wide cells against the masses of the
// truncated Gaussian, and the number of cells used.
double spectrumChiSquare(const std::vector<double>& En, double mu, double sigma,
                         std::size_t& cells) {
    const double total = st::normalMassOnDomain(kDomain, mu, sigma);
    const double n = static_cast<double>(En.size());
    double chi2 = 0.0;
    cells = 0;
    for (double lo = kDomain.En_min; lo < kDomain.En_max - 1e-9; lo += 1.0) {
        const double hi = std::min(lo + 1.0, kDomain.En_max);
        const double expected =
            n * st::normalMassOnDomain({lo, hi, st::kEnergyGridPoints}, mu, sigma) / total;
        const auto observed = static_cast<double>(
            std::count_if(En.begin(), En.end(), [&](double e) { return e >= lo && e < hi; }));
        if (expected > 5.0) {
            chi2 += (observed - expected) * (observed - expected) / expected;
            ++cells;
        }
    }
    return chi2;
}

} // namespace

TEST_CASE("the pool histogram is a density on the domain", "[resample_c]") {
    const std::vector<double> pool = risingPool(200000, kSeed);
    const st::BinnedDensity g = st::histogramDensity(pool, kDomain);
    REQUIRE(g.density.size() == 95);
    double integral = 0.0;
    for (const double d : g.density) {
        integral += d * g.step;
    }
    CHECK_THAT(integral, WithinRel(1.0, 1e-12));
    CHECK(g.binOf(kDomain.En_min) == 0);
    CHECK(g.binOf(kDomain.En_max) == 94);
    CHECK(g.binOf(10.05) == 0);
    CHECK_THROWS_AS((void)g.binOf(9.99), std::invalid_argument);
    CHECK_THROWS_AS((void)g.binOf(19.51), std::invalid_argument);
}

TEST_CASE("the pool histogram rejects what would break the rejection", "[resample_c]") {
    // An empty bin: no events above 19 MeV.
    std::vector<double> pool = risingPool(10000, kSeed);
    for (double& En : pool) {
        En = std::min(En, 18.99);
    }
    CHECK_THROWS_AS((void)st::histogramDensity(pool, kDomain), std::invalid_argument);
    CHECK_THROWS_AS((void)st::histogramDensity({9.0}, kDomain), std::invalid_argument);
    // 0.3 MeV does not tile 9.5 MeV.
    CHECK_THROWS_AS((void)st::histogramDensity(risingPool(10000, kSeed), kDomain, 0.3),
                    std::invalid_argument);
}

TEST_CASE("the truncated density integrates to one and vanishes outside", "[resample_c]") {
    constexpr std::size_t kSteps = 100000;
    const double step = (kDomain.En_max - kDomain.En_min) / static_cast<double>(kSteps);
    double integral = 0.0;
    for (std::size_t s = 0; s < kSteps; ++s) {
        const double En = kDomain.En_min + ((static_cast<double>(s) + 0.5) * step);
        integral += st::truncatedNormalDensity(kDomain, 18.0, 2.4, En) * step;
    }
    CHECK_THAT(integral, WithinRel(1.0, 1e-6));
    CHECK_THAT(st::truncatedNormalDensity(kDomain, 15.0, 1.0, 9.0), WithinAbs(0.0, 0.0));
    CHECK_THROWS_AS((void)st::truncatedNormalDensity(kDomain, 40.0, 1.0, 15.0),
                    std::invalid_argument);
}

TEST_CASE("the resampling bound is the maximum of f over g_hat", "[resample_c]") {
    const st::BinnedDensity g = st::histogramDensity(risingPool(200000, kSeed), kDomain);
    for (const auto& [mu, sigma] : {std::pair{13.0, 0.8}, std::pair{17.0, 2.4},
                                    std::pair{15.03, 1.1}, std::pair{19.4, 0.8}}) {
        const double c = st::resampleBound(g, mu, sigma);
        // Both densities integrate to 1, so f <= c g needs c >= 1.
        CHECK(c >= 1.0);
        // A dense scan never exceeds c and comes close to it. The maximum may
        // sit on the upper edge of a bin, which belongs to the next bin, so
        // the scan reaches it only to within one step (1e-5 MeV); there
        // d ln f / dE = (mu - E) / sigma^2 is at most ~10 / MeV, hence 1e-4.
        double scan = 0.0;
        constexpr std::size_t kScan = 950001;
        for (std::size_t s = 0; s < kScan; ++s) {
            const double En =
                kDomain.En_min + ((kDomain.En_max - kDomain.En_min) * static_cast<double>(s) /
                                  static_cast<double>(kScan - 1));
            scan = std::max(scan, st::truncatedNormalDensity(kDomain, mu, sigma, En) /
                                      g.density[g.binOf(En)]);
        }
        INFO("mu = " << mu << " sigma = " << sigma);
        CHECK(scan <= c * (1.0 + 1e-12));
        CHECK_THAT(scan, WithinRel(c, 1e-4));
    }
}

TEST_CASE("the bound against a flat pool is f at its peak over the flat density", "[resample_c]") {
    // One event per bin: g_hat = 1 / 9.5 MeV everywhere, so c = f(mu) * 9.5.
    // A very wide f is flat on the domain, where c -> 1 from above.
    std::vector<double> flat;
    for (std::size_t b = 0; b < 95; ++b) {
        flat.push_back(kDomain.En_min + (0.1 * (static_cast<double>(b) + 0.5)));
    }
    const st::BinnedDensity g = st::histogramDensity(flat, kDomain);
    const double wide = st::resampleBound(g, 14.75, 1e3);
    CHECK(wide >= 1.0);
    CHECK_THAT(wide, WithinRel(1.0, 1e-4));
    const double expected = st::truncatedNormalDensity(kDomain, 14.75, 1.0, 14.75) * 9.5;
    CHECK_THAT(st::resampleBound(g, 14.75, 1.0), WithinRel(expected, 1e-12));
}

TEST_CASE("resampled experiments follow f and are disjoint", "[resample_c]") {
    const std::vector<double> pool = risingPool(400000, kSeed);
    const st::BinnedDensity g = st::histogramDensity(pool, kDomain);
    constexpr std::size_t kN = 1000;
    constexpr std::size_t kM = 40;
    for (const auto& [mu, sigma] : {std::pair{13.0, 0.8}, std::pair{16.5, 2.4}}) {
        auto engine = rnd::perEventEngine(kSeed, 1);
        const st::ResampledRow row =
            st::resampleRow(pool, g, sameTruths(kM, mu, sigma), kN, engine);
        REQUIRE_FALSE(row.exhausted());
        REQUIRE(row.experiments.size() == kM);
        REQUIRE(row.visited.size() == kM);
        REQUIRE(row.bounds.size() == kM);
        std::set<std::size_t> seen;
        std::vector<double> kept;
        std::size_t visited = 0;
        for (std::size_t j = 0; j < kM; ++j) {
            REQUIRE(row.experiments[j].size() == kN);
            visited += row.visited[j];
            for (const std::size_t index : row.experiments[j]) {
                REQUIRE(seen.insert(index).second);
                kept.push_back(pool[index]);
            }
        }
        CHECK(visited == row.visited_total);
        // The mean number of visits per kept event is c.
        const double ratio = static_cast<double>(visited) / static_cast<double>(kM * kN);
        CHECK_THAT(ratio, WithinRel(row.bounds.front(), 0.05));
        std::size_t cells = 0;
        const double chi2 = spectrumChiSquare(kept, mu, sigma, cells);
        INFO("mu = " << mu << " chi2 = " << chi2 << " on " << cells << " cells");
        REQUIRE(cells >= 4);
        // chi^2 well below 3 cells' worth above the dof: P < 1e-3 for <= 10
        // cells. The seed is fixed, so the check is deterministic.
        CHECK(chi2 < static_cast<double>(cells) + 25.0);
        // The pool itself is clearly rejected: the test has power.
        std::vector<double> raw(pool.begin(),
                                pool.begin() + static_cast<std::ptrdiff_t>(kept.size()));
        std::size_t raw_cells = 0;
        CHECK(spectrumChiSquare(raw, mu, sigma, raw_cells) > 1e3);
    }
}

TEST_CASE("a short pool keeps only the complete experiments", "[resample_c]") {
    const std::vector<double> pool = risingPool(20000, kSeed);
    const st::BinnedDensity g = st::histogramDensity(pool, kDomain);
    auto engine = rnd::perEventEngine(kSeed, 2);
    const st::ResampledRow row = st::resampleRow(pool, g, sameTruths(40, 15.0, 1.0), 1000, engine);
    CHECK(row.planned == 40);
    CHECK(row.exhausted());
    REQUIRE_FALSE(row.experiments.empty());
    CHECK(row.experiments.size() < 40);
    CHECK(row.visited.size() == row.experiments.size());
    for (const auto& experiment : row.experiments) {
        CHECK(experiment.size() == 1000);
    }
    CHECK(row.visited_total == pool.size());
    CHECK_THROWS_AS((void)st::resampleRow(pool, g, sameTruths(1, 15.0, 1.0), 0, engine),
                    std::invalid_argument);
}

TEST_CASE("a resampled row depends only on the pool, truths and engine", "[resample_c]") {
    const std::vector<double> pool = risingPool(50000, kSeed);
    const st::BinnedDensity g = st::histogramDensity(pool, kDomain);
    auto e1 = rnd::perEventEngine(kSeed, 3);
    auto e2 = rnd::perEventEngine(kSeed, 3);
    auto e3 = rnd::perEventEngine(kSeed, 4);
    const auto truths = sameTruths(5, 14.0, 1.5);
    const st::ResampledRow a = st::resampleRow(pool, g, truths, 300, e1);
    const st::ResampledRow b = st::resampleRow(pool, g, truths, 300, e2);
    const st::ResampledRow c = st::resampleRow(pool, g, truths, 300, e3);
    CHECK(a.experiments == b.experiments);
    CHECK(a.visited == b.visited);
    CHECK(a.experiments != c.experiments);
}

TEST_CASE("the truth at a fixed direction is the one of the toy experiment", "[resample_c]") {
    const geo::Direction omega{1.0, 0.0, 0.0};
    for (const std::uint64_t index : {0U, 7U, 199U}) {
        const st::CaseCPseudoExperiment toy =
            st::caseCPseudoExperiment(kDomain, 50, kResolution, st::kCaseCChecklistTruth,
                                      st::CaseCSpectrum::truncated, kSeed, 16, index);
        const st::CaseCTruth truth =
            st::caseCTruthAtDirection(st::kCaseCChecklistTruth, omega, kSeed, 16, index);
        CHECK_THAT(truth.mu, WithinULP(toy.truth.mu, 0));
        CHECK_THAT(truth.sigma, WithinULP(toy.truth.sigma, 0));
        CHECK_THAT(truth.omega_n.x, WithinULP(1.0, 0));
        // At the toy's own truth the twin draws the very same events.
        const st::CaseCPseudoExperiment twin =
            st::caseCPseudoExperimentAtTruth(kDomain, 50, kResolution, toy.truth, kSeed, 16, index);
        for (std::size_t k = 0; k < 50; ++k) {
            REQUIRE_THAT(twin.sample.truth_En[k], WithinULP(toy.sample.truth_En[k], 0));
            REQUIRE_THAT(twin.sample.reco_Ep[k], WithinULP(toy.sample.reco_Ep[k], 0));
            REQUIRE_THAT(twin.sample.reco_track[k].z, WithinULP(toy.sample.reco_track[k].z, 0));
        }
    }
}

TEST_CASE("the accepted toy twin keeps the E_n of the ideal one", "[resample_c]") {
    st::AcceptanceHistogram h;
    for (std::size_t k = 0; k <= 10; ++k) {
        h.En_edges.push_back(10.0 + static_cast<double>(k));
    }
    constexpr std::size_t kThetaBins = 9;
    for (std::size_t t = 0; t <= kThetaBins; ++t) {
        h.theta_edges.push_back(0.5 * kPi * static_cast<double>(t) /
                                static_cast<double>(kThetaBins));
    }
    h.theta_edges.back() = 0.5 * kPi;
    for (std::size_t k = 0; k < 10; ++k) {
        for (std::size_t t = 0; t < kThetaBins; ++t) {
            // A forward deficit, so the accept loop really rejects.
            h.counts.push_back(t < 3 ? 0.5 : 1.0);
        }
    }
    const st::AcceptanceRatio ratio(h);
    const st::CaseCTruth truth =
        st::caseCTruthAtDirection(st::kCaseCChecklistTruth, kOmega, kSeed, 16, 3);
    const auto ideal =
        st::caseCPseudoExperimentAtTruth(kDomain, 500, kResolution, truth, kSeed, 16, 3);
    const auto accepted =
        st::caseCPseudoExperimentAtTruth(kDomain, 500, kResolution, truth, kSeed, 16, 3, &ratio);
    std::size_t differ = 0;
    for (std::size_t k = 0; k < 500; ++k) {
        REQUIRE_THAT(accepted.sample.truth_En[k], WithinULP(ideal.sample.truth_En[k], 0));
        REQUIRE_THAT(geo::angleBetween(accepted.sample.truth_track[k], kOmega),
                     WithinAbs(accepted.sample.truth_theta_p[k], 1e-9));
        if (std::abs(accepted.sample.truth_theta_p[k] - ideal.sample.truth_theta_p[k]) > 0.0) {
            ++differ;
        }
    }
    // Some events were redrawn: the forward bins keep half their tries.
    CHECK(differ > 0);
}

TEST_CASE("the Case C model passes its acceptance to stage 1", "[resample_c]") {
    const st::Grid En = st::energyGrid({});
    const st::Grid theta = st::recoilAngleGrid(st::kTrackAngleGridPoints);
    const st::Grid theta_obs{0.0, kPi, st::kTrackObservedGridPoints};
    const st::CaseCModel ideal(En, theta, theta_obs, kResolution);
    // A constant log acceptance shifts every energy posterior by that
    // constant: the base moves by it and nothing else.
    constexpr double kShift = 0.3;
    const std::vector<double> log_acceptance(En.size() * theta.size(), kShift);
    const st::CaseCModel shifted(En, theta, theta_obs, kResolution, log_acceptance);
    CHECK_FALSE(ideal.stage1Model().hasAcceptance());
    REQUIRE(shifted.stage1Model().hasAcceptance());
    CHECK(shifted.stage1Model().logAcceptance().size() == log_acceptance.size());
    const std::vector<double> Ep_hat{11.2, 6.8, 13.9};
    std::vector<geo::Direction> tracks;
    for (const double t : {0.52, 0.86, 0.27}) {
        tracks.push_back({std::sin(t), 0.0, std::cos(t)});
    }
    const st::HierarchicalBase a = ideal.base(Ep_hat, tracks, kOmega);
    const st::HierarchicalBase b = shifted.base(Ep_hat, tracks, kOmega);
    REQUIRE(a.values.size() == b.values.size());
    for (std::size_t i = 0; i < a.values.size(); ++i) {
        if (std::isfinite(a.values[i])) {
            REQUIRE_THAT(b.values[i] - a.values[i], WithinAbs(kShift, 1e-9));
        }
    }
}
