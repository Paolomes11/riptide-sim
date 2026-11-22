#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "riptide/common/random.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/validate.hpp"

namespace st = riptide::stats;
namespace rnd = riptide::random;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr std::uint64_t kSeed = 20260907;

double normalCdf(double z) {
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

// Standard normal quantile by bisection on the CDF, so the test does not
// hard-code it.
double normalQuantile(double p) {
    double lo = -10.0;
    double hi = 10.0;
    constexpr int kIterations = 200;
    for (int i = 0; i < kIterations; ++i) {
        const double mid = 0.5 * (lo + hi);
        (normalCdf(mid) < p ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

std::vector<double> gaussianLogPost(const st::Grid& grid, double mean, double sigma) {
    std::vector<double> lp(grid.size());
    for (std::size_t i = 0; i < grid.size(); ++i) {
        const double u = (grid[i] - mean) / sigma;
        lp[i] = -0.5 * u * u;
    }
    return lp;
}

} // namespace

TEST_CASE("posterior moments of a Gaussian on a fine grid", "[validate]") {
    // Symmetric about the mean, so the truncation does not shift it.
    const st::Grid grid(-13.0, 17.0, 3001);
    const auto lp = gaussianLogPost(grid, 2.0, 1.5);
    const auto m = st::posteriorMeanStd(grid, lp);
    CHECK_THAT(m.mean, WithinAbs(2.0, 1e-9));
    CHECK_THAT(m.std, WithinRel(1.5, 1e-4));
    // An additive constant does not change normalised weights.
    auto shifted = lp;
    for (double& v : shifted) {
        v += 700.0;
    }
    const auto ms = st::posteriorMeanStd(grid, shifted);
    CHECK_THAT(ms.mean, WithinAbs(m.mean, 1e-12));
    CHECK_THAT(ms.std, WithinAbs(m.std, 1e-12));
}

TEST_CASE("credible intervals of a Gaussian are within one grid step", "[validate]") {
    const st::Grid grid(-10.0, 10.0, 401);
    const double mean = 1.3;
    const double sigma = 1.7;
    const auto lp = gaussianLogPost(grid, mean, sigma);
    const std::vector<double> levels{0.68, 0.90, 0.95};
    const auto intervals = st::credibleIntervals(grid, lp.data(), levels);
    REQUIRE(intervals.size() == levels.size());
    for (std::size_t l = 0; l < levels.size(); ++l) {
        const double z = normalQuantile(0.5 * (1.0 + levels[l]));
        CHECK_THAT(intervals[l].lo, WithinAbs(mean - (z * sigma), grid.step()));
        CHECK_THAT(intervals[l].hi, WithinAbs(mean + (z * sigma), grid.step()));
    }
    CHECK(intervals[0].lo > intervals[1].lo);
    CHECK(intervals[1].lo > intervals[2].lo);
}

TEST_CASE("credible intervals stop at the grid edges", "[validate]") {
    const st::Grid grid(0.0, 1.0, 11);
    std::vector<double> lp(grid.size(), -1000.0);
    lp[0] = 0.0; // all the mass on the first point
    const auto iv = st::credibleIntervals(grid, lp.data(), {0.9});
    CHECK(iv[0].lo >= grid.front());
    CHECK(iv[0].hi <= grid[1]);
}

TEST_CASE("validation inputs are checked", "[validate]") {
    const st::Grid grid(0.0, 1.0, 5);
    const std::vector<double> lp(grid.size(), 0.0);
    CHECK_THROWS_AS(st::credibleIntervals(grid, lp.data(), {1.0}), std::invalid_argument);
    CHECK_THROWS_AS(st::credibleIntervals(grid, lp.data(), {0.0}), std::invalid_argument);
    const std::vector<double> all_minus_inf(grid.size(), -INFINITY);
    CHECK_THROWS_AS(st::posteriorMeanStd(grid, all_minus_inf), std::invalid_argument);
    CHECK_THROWS_AS(st::posteriorMeanStd(grid, std::vector<double>(3, 0.0)), std::invalid_argument);
    CHECK_THROWS_AS(st::pulls({1.0}, {1.0}, {0.0}), std::invalid_argument);
    CHECK_THROWS_AS(st::pulls({1.0}, {1.0, 2.0}, {1.0}), std::invalid_argument);
    CHECK_THROWS_AS(st::biasCurve({}, {}, 3), std::invalid_argument);
    CHECK_THROWS_AS(st::biasCurve({1.0, 2.0}, {1.0, 2.0}, 0), std::invalid_argument);
}

TEST_CASE("coverage includes the interval end points", "[validate]") {
    st::IntervalTable table;
    table.levels = {0.5};
    table.values = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}};
    const auto c = st::coverage({1.0, 2.0, 1.5, 2.5}, table);
    REQUIRE(c.size() == 1);
    CHECK_THAT(c[0], WithinAbs(0.75, 1e-15));
}

TEST_CASE("bias curve binning follows numpy.digitize", "[validate]") {
    // Edges 0, 1, 2, 3: 1 goes to bin 1, the maximum 3 to the last bin and
    // bin 2 has an event only through the maximum.
    const std::vector<double> truth{0.0, 0.5, 1.0, 3.0};
    const std::vector<double> estimate{0.1, 0.7, 1.4, 3.8};
    const auto curve = st::biasCurve(truth, estimate, 3);
    REQUIRE(curve.counts == std::vector<std::size_t>{2, 1, 1});
    CHECK_THAT(curve.centers[0], WithinAbs(0.5, 1e-15));
    CHECK_THAT(curve.bias[0], WithinAbs(0.15, 1e-12));
    CHECK_THAT(curve.bias[1], WithinAbs(0.4, 1e-12));
    CHECK_THAT(curve.bias[2], WithinAbs(0.8, 1e-12));

    const auto sparse = st::biasCurve({0.0, 3.0}, {0.0, 3.0}, 3);
    CHECK(std::isnan(sparse.bias[1]));
    CHECK(sparse.counts[1] == 0);
}

TEST_CASE("R3: the checklist sees a known miscalibration", "[validate][oracle]") {
    // Estimator with a multiplicative and additive bias and a declared sigma
    // 30% too small: the checklist must report each defect.
    constexpr std::size_t kEvents = 20000;
    constexpr double kTruthLo = 1.0;
    constexpr double kTruthHi = 5.0;
    constexpr double kRelSigma = 0.08;
    constexpr double kScale = 1.03;
    constexpr double kOffset = 0.05;
    constexpr double kSigmaUnderstatement = 0.7;
    constexpr double kLevel = 0.68;

    std::vector<double> truth(kEvents);
    std::vector<double> estimate(kEvents);
    std::vector<double> sigma_hat(kEvents);
    st::IntervalTable intervals;
    intervals.levels = {kLevel};
    intervals.values.resize(kEvents);
    const double z = normalQuantile(0.5 * (1.0 + kLevel));
    for (std::size_t k = 0; k < kEvents; ++k) {
        auto engine = rnd::perEventEngine(kSeed, k);
        truth[k] = kTruthLo + ((kTruthHi - kTruthLo) * rnd::uniform53(engine));
        const double true_sigma = kRelSigma * truth[k];
        estimate[k] = (kScale * truth[k]) + kOffset + (true_sigma * rnd::standardNormal(engine));
        sigma_hat[k] = kSigmaUnderstatement * true_sigma;
        intervals.values[k] = {estimate[k] - (z * sigma_hat[k]), estimate[k] + (z * sigma_hat[k])};
    }

    const auto result = st::runChecklist(truth, estimate, sigma_hat, intervals, 20);
    CHECK_THAT(result.bias.bias.front(), WithinAbs(0.08, 0.02));
    CHECK_THAT(result.bias.bias.back(), WithinAbs(0.20, 0.02));
    CHECK_THAT(result.pull_mean, WithinAbs(0.9, 0.1));
    CHECK_THAT(result.pull_width, WithinAbs(1.4, 0.1));
    REQUIRE(result.coverage.size() == 1);
    CHECK_THAT(result.coverage[0], WithinAbs(0.45, 0.05));
}
