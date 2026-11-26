#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "riptide/stats/case_c.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/log_sum_exp.hpp"

namespace st = riptide::stats;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Trapezoidal integral of exp(log_density) on the grid.
double trapezoidIntegral(const st::Grid& grid, const std::vector<double>& log_density) {
    double sum = 0.0;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        const double weight = (i == 0 || i + 1 == grid.size()) ? 0.5 : 1.0;
        sum += weight * std::exp(log_density[i]);
    }
    return sum * grid.step();
}

double logMassSum(const std::vector<double>& log_mass) {
    return std::exp(st::logSumExp(log_mass));
}

} // namespace

TEST_CASE("the truncated Gaussian prior integrates to 1 on the domain", "[case_c]") {
    const st::Grid En = st::energyGrid({});
    constexpr double kTolerance = 1e-12;
    // mu inside, on an edge and outside the domain; sigma from below the
    // grid step to far above the domain width.
    for (const double mu : {15.0, 10.0, 20.0, 8.0, 23.0}) {
        for (const double sigma : {0.01, 1.0, 50.0}) {
            const std::vector<double> log_prior = st::truncatedNormalEnergyLogPrior(En, mu, sigma);
            REQUIRE(log_prior.size() == En.size());
            CHECK_THAT(trapezoidIntegral(En, log_prior), WithinAbs(1.0, kTolerance));
        }
    }
}

TEST_CASE("a very wide truncated Gaussian is the flat prior", "[case_c]") {
    const st::Grid En = st::energyGrid({});
    // At sigma = 1e6 MeV the Gaussian varies by (10 / 1e6)^2 / 2 = 5e-11 on
    // the domain.
    const std::vector<double> log_prior = st::truncatedNormalEnergyLogPrior(En, 15.0, 1e6);
    const double flat = -std::log(En.back() - En.front());
    for (const double value : log_prior) {
        CHECK_THAT(value, WithinRel(flat, 1e-9));
    }
}

TEST_CASE("a narrow truncated Gaussian stays finite far from its mean", "[case_c]") {
    const st::Grid En = st::energyGrid({});
    const std::vector<double> log_prior = st::truncatedNormalEnergyLogPrior(En, 12.0, 0.01);
    for (const double value : log_prior) {
        CHECK(std::isfinite(value));
    }
    // 8 MeV from the mean, i.e. 800 sigma: exp would underflow, the log not.
    CHECK(log_prior.back() < -3e5);
}

TEST_CASE("the truncated Gaussian prior rejects invalid parameters", "[case_c]") {
    const st::Grid En = st::energyGrid({});
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(st::truncatedNormalEnergyLogPrior(En, 15.0, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(st::truncatedNormalEnergyLogPrior(En, 15.0, -1.0), std::invalid_argument);
    CHECK_THROWS_AS(st::truncatedNormalEnergyLogPrior(En, nan, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(st::truncatedNormalEnergyLogPrior(En, 15.0, nan), std::invalid_argument);
}

TEST_CASE("the global hyperparameter grid spans the domain and the sigma range", "[case_c]") {
    const st::HyperGrid grid = st::globalHyperGrid({});
    REQUIRE(grid.size() == st::kHyperMuPoints * st::kHyperSigmaPoints);
    CHECK_THAT(grid.mu.front(), WithinRel(10.0, 1e-15));
    CHECK_THAT(grid.mu.back(), WithinRel(20.0, 1e-15));
    CHECK_THAT(grid.sigmaAt(0), WithinRel(st::kHyperSigmaMin, 1e-14));
    CHECK_THAT(grid.sigmaAt(grid.size() - 1), WithinRel(st::kHyperSigmaMax, 1e-14));
    // mu outer, sigma inner.
    const std::size_t index = grid.index(3, 7);
    CHECK_THAT(grid.muAt(index), WithinRel(grid.mu[3], 1e-15));
    CHECK_THAT(grid.sigmaAt(index), WithinRel(std::exp(grid.log_sigma[7]), 1e-15));
    // Geometric sigma: a constant ratio between neighbours.
    const double ratio = grid.sigmaAt(1) / grid.sigmaAt(0);
    CHECK_THAT(grid.sigmaAt(31) / grid.sigmaAt(30), WithinRel(ratio, 1e-12));
}

TEST_CASE("the window holds every point near the maximum, one step wider", "[case_c]") {
    const st::HyperGrid grid = st::globalHyperGrid({});
    // A Gaussian bump in (mu, log sigma) centred inside the grid.
    constexpr double kMu0 = 14.3;
    constexpr double kLogSigma0 = 0.2;
    constexpr double kWidthMu = 0.4;
    constexpr double kWidthLogSigma = 0.3;
    std::vector<double> log_post(grid.size());
    for (std::size_t i = 0; i < grid.mu.size(); ++i) {
        for (std::size_t j = 0; j < grid.log_sigma.size(); ++j) {
            const double a = (grid.mu[i] - kMu0) / kWidthMu;
            const double b = (grid.log_sigma[j] - kLogSigma0) / kWidthLogSigma;
            log_post[grid.index(i, j)] = -0.5 * ((a * a) + (b * b));
        }
    }
    const st::HyperWindow window = st::hyperWindow(log_post, grid);
    const double max_value = 0.0; // bounded above by the bump peak
    bool max_inside = false;
    double best = -std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < grid.size(); ++index) {
        const double mu = grid.muAt(index);
        const double log_sigma = std::log(grid.sigmaAt(index));
        const bool inside = mu > window.mu_lo && mu < window.mu_hi &&
                            log_sigma > window.log_sigma_lo && log_sigma < window.log_sigma_hi;
        if (log_post[index] > max_value - st::kHyperWindowDeltaLog) {
            // Strictly inside: the window is one step wider than the kept points.
            CHECK(inside);
        }
        if (log_post[index] > best) {
            best = log_post[index];
            max_inside = inside;
        }
    }
    CHECK(max_inside);
    const st::HyperGrid fine = st::windowHyperGrid(window);
    CHECK(fine.size() == st::kHyperFinePoints * st::kHyperFinePoints);
    CHECK_THAT(fine.mu.front(), WithinRel(window.mu_lo, 1e-15));
    CHECK_THAT(fine.log_sigma.back(), WithinRel(window.log_sigma_hi, 1e-15));
}

TEST_CASE("the window is clipped to the global grid", "[case_c]") {
    const st::HyperGrid grid = st::globalHyperGrid({});
    // The maximum on the lowest corner of the grid.
    std::vector<double> log_post(grid.size());
    for (std::size_t i = 0; i < grid.mu.size(); ++i) {
        for (std::size_t j = 0; j < grid.log_sigma.size(); ++j) {
            log_post[grid.index(i, j)] = -static_cast<double>(i + j);
        }
    }
    const st::HyperWindow window = st::hyperWindow(log_post, grid);
    CHECK_THAT(window.mu_lo, WithinRel(grid.mu.front(), 1e-15));
    CHECK_THAT(window.log_sigma_lo, WithinRel(grid.log_sigma.front(), 1e-15));
    // Points with i + j < 10 are kept: index 9 at most, plus one step.
    CHECK_THAT(window.mu_hi, WithinRel(grid.mu[10], 1e-12));
    CHECK_THAT(window.log_sigma_hi, WithinRel(grid.log_sigma[10], 1e-12));
}

TEST_CASE("the window rejects a log-posterior without a finite maximum", "[case_c]") {
    const st::HyperGrid grid = st::globalHyperGrid({});
    const std::vector<double> minus_inf(grid.size(), -std::numeric_limits<double>::infinity());
    CHECK_THROWS_AS(st::hyperWindow(minus_inf, grid), std::invalid_argument);
    CHECK_THROWS_AS(st::hyperWindow(std::vector<double>(3, 0.0), grid), std::invalid_argument);
}

TEST_CASE("the hyperparameter priors are proper", "[case_c]") {
    const st::HyperGrid grid = st::globalHyperGrid({});
    const std::vector<double> flat = st::uniformLogSigmaHyperPrior(grid);
    const std::vector<double> linear = st::uniformSigmaHyperPrior(grid);
    CHECK_THAT(logMassSum(flat), WithinAbs(1.0, 1e-12));
    CHECK_THAT(logMassSum(linear), WithinAbs(1.0, 1e-12));
    // Uniform in sigma: the mass ratio between two points is the sigma ratio.
    const std::size_t a = grid.index(5, 2);
    const std::size_t b = grid.index(40, 50);
    CHECK_THAT(std::exp(linear[b] - linear[a]),
               WithinRel(grid.sigmaAt(b) / grid.sigmaAt(a), 1e-12));
    CHECK_THAT(flat[a], WithinRel(flat[b], 1e-15));
}
