#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/sim/check_stats.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace checks = riptide::sim::checks;

namespace {

// Common seed of the project, so that the sampled tests are reproducible.
constexpr std::uint64_t kSeed = 20260907;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// Normalized f(mu) = sum_l (2l+1)/2 a_l P_l(mu), a_0 = 1.
double legendreDensity(const std::vector<double>& a, double mu) {
    double value = 0.5;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const int degree = static_cast<int>(i) + 1;
        value += 0.5 * (2.0 * degree + 1.0) * a[i] * riptide::kinematics::legendreP(degree, mu);
    }
    return value;
}

// Upper bound of f on [-1, 1], since |P_l| <= 1 there.
double densityBound(const std::vector<double>& a) {
    double bound = 0.5;
    for (std::size_t i = 0; i < a.size(); ++i) {
        bound += 0.5 * (2.0 * static_cast<double>(i + 1) + 1.0) * std::abs(a[i]);
    }
    return bound;
}

std::vector<double> sampleLegendre(const std::vector<double>& a, std::size_t count,
                                   std::mt19937_64& rng) {
    std::uniform_real_distribution<double> cosine(-1.0, 1.0);
    std::uniform_real_distribution<double> height(0.0, densityBound(a));
    std::vector<double> mu;
    mu.reserve(count);
    while (mu.size() < count) {
        const double candidate = cosine(rng);
        if (height(rng) < legendreDensity(a, candidate)) {
            mu.push_back(candidate);
        }
    }
    return mu;
}

// Forward-peaked but positive everywhere, like the n-p distribution at 20 MeV
// only much more anisotropic, so that every coefficient is well measured.
const std::vector<double> kCoefficients = {0.3, 0.1, -0.05, 0.02};

} // namespace

TEST_CASE("Pearson chi2 of a histogram", "[checks][chi2]") {
    const checks::HistogramComparison result =
        checks::compareHistogram({10, 20, 30}, {12, 18, 30}, 1);
    CHECK(result.ndf == 2);
    CHECK_THAT(result.chi2, WithinRel(4.0 / 12.0 + 4.0 / 18.0, 1e-14));
    REQUIRE(result.pulls.size() == 3);
    CHECK_THAT(result.pulls[0], WithinRel(-2.0 / std::sqrt(12.0), 1e-14));
    CHECK_THAT(result.pulls[1], WithinRel(2.0 / std::sqrt(18.0), 1e-14));
    CHECK_THAT(result.pulls[2], WithinAbs(0.0, 0.0));
    CHECK_THAT(result.min_expected, WithinAbs(12.0, 0.0));
    CHECK_THAT(result.p_value, WithinRel(std::exp(-0.5 * result.chi2), 1e-13));

    CHECK_THROWS_AS(checks::compareHistogram({1, 2}, {6, 6, 6}, 0), std::invalid_argument);
    CHECK_THROWS_AS(checks::compareHistogram({10, 10}, {10, 4.9}, 0), std::invalid_argument);
    CHECK_THROWS_AS(checks::compareHistogram({-1, 10}, {10, 10}, 0), std::invalid_argument);
    CHECK_THROWS_AS(checks::compareHistogram({10, 10}, {10, 10}, 2), std::invalid_argument);
    CHECK_THROWS_AS(checks::compareHistogram({kNaN, 10}, {10, 10}, 0), std::invalid_argument);
}

TEST_CASE("chi2 survival probability against closed forms", "[checks][chi2]") {
    // The values cover both the series (x < a + 1) and the continued fraction.
    for (const double chi2 : {0.01, 0.5, 1.0, 2.0, 3.9, 4.1, 8.0, 30.0, 200.0}) {
        CAPTURE(chi2);
        // ndf = 2: Q = e^{-chi2/2}.
        CHECK_THAT(checks::chi2SurvivalProbability(chi2, 2),
                   WithinRel(std::exp(-0.5 * chi2), 1e-12));
        // ndf = 4: Q = e^{-chi2/2} (1 + chi2/2).
        CHECK_THAT(checks::chi2SurvivalProbability(chi2, 4),
                   WithinRel(std::exp(-0.5 * chi2) * (1.0 + 0.5 * chi2), 1e-12));
        // ndf = 1: Q = erfc(sqrt(chi2/2)).
        CHECK_THAT(checks::chi2SurvivalProbability(chi2, 1),
                   WithinRel(std::erfc(std::sqrt(0.5 * chi2)), 1e-12));
    }
    CHECK_THAT(checks::chi2SurvivalProbability(0.0, 3), WithinAbs(1.0, 0.0));
    CHECK_THROWS_AS(checks::chi2SurvivalProbability(1.0, 0), std::invalid_argument);
    CHECK_THROWS_AS(checks::chi2SurvivalProbability(-1.0, 3), std::invalid_argument);
    CHECK_THROWS_AS(checks::chi2SurvivalProbability(kInf, 3), std::invalid_argument);
}

TEST_CASE("chi2 survival probability is continuous at the series/fraction switch",
          "[checks][chi2]") {
    for (const int ndf : {1, 3, 10, 51, 400}) {
        CAPTURE(ndf);
        // The switch is at chi2/2 = ndf/2 + 1.
        const double edge = ndf + 2.0;
        const double below = checks::chi2SurvivalProbability(std::nextafter(edge, 0.0), ndf);
        const double above = checks::chi2SurvivalProbability(edge, ndf);
        CHECK_THAT(below, WithinRel(above, 1e-12));
    }
}

TEST_CASE("fractions of categories against their probabilities", "[checks][fractions]") {
    SECTION("exact agreement") {
        const checks::FractionComparison result =
            checks::compareFractions({50, 30, 20}, {0.5, 0.3, 0.2});
        CHECK(result.total == 100);
        CHECK(result.ndf == 2);
        CHECK_THAT(result.chi2, WithinAbs(0.0, 1e-24));
        CHECK_THAT(result.p_value, WithinAbs(1.0, 1e-15));
        for (const double pull : result.pulls) {
            CHECK_THAT(pull, WithinAbs(0.0, 1e-12));
        }
    }
    SECTION("hand-computed discrepancy") {
        const checks::FractionComparison result =
            checks::compareFractions({60, 25, 15}, {0.5, 0.3, 0.2});
        CHECK_THAT(result.chi2, WithinRel(100.0 / 50.0 + 25.0 / 30.0 + 25.0 / 20.0, 1e-14));
        CHECK_THAT(result.fractions[0], WithinRel(0.6, 1e-15));
        CHECK_THAT(result.sigmas[0], WithinRel(0.05, 1e-14));
        CHECK_THAT(result.pulls[0], WithinRel(2.0, 1e-13));
    }
    SECTION("an empty category keeps a finite sigma") {
        const checks::FractionComparison result = checks::compareFractions({0, 100}, {0.05, 0.95});
        CHECK_THAT(result.fractions[0], WithinAbs(0.0, 0.0));
        CHECK_THAT(result.pulls[0], WithinRel(-0.05 / std::sqrt(0.05 * 0.95 / 100.0), 1e-13));
    }
    SECTION("invalid input") {
        CHECK_THROWS_AS(checks::compareFractions({10}, {1.0}), std::invalid_argument);
        CHECK_THROWS_AS(checks::compareFractions({10, 10}, {0.5, 0.4}), std::invalid_argument);
        CHECK_THROWS_AS(checks::compareFractions({10, 10}, {0.0, 1.0}), std::invalid_argument);
        CHECK_THROWS_AS(checks::compareFractions({50, 50}, {0.04, 0.96}), std::invalid_argument);
        CHECK_THROWS_AS(checks::compareFractions({50, 50}, {0.5, 0.5, 0.0}), std::invalid_argument);
    }
}

TEST_CASE("Legendre moments recover the coefficients of a sample", "[checks][legendre]") {
    // The test density must be a density: positive on the whole interval.
    constexpr int kGridPoints = 2001;
    for (int i = 0; i < kGridPoints; ++i) {
        const double mu = -1.0 + 2.0 * i / (kGridPoints - 1.0);
        REQUIRE(legendreDensity(kCoefficients, mu) > 0.0);
    }
    std::mt19937_64 rng(kSeed);
    constexpr int kDegrees = 6; // two beyond the sampled ones, which must vanish

    SECTION("one large sample") {
        const checks::LegendreMoments moments =
            checks::legendreMoments(sampleLegendre(kCoefficients, 200000, rng), kDegrees);
        CHECK(moments.samples == 200000);
        REQUIRE(moments.values.size() == kDegrees);
        for (std::size_t i = 0; i < moments.values.size(); ++i) {
            CAPTURE(i + 1);
            const double truth = i < kCoefficients.size() ? kCoefficients[i] : 0.0;
            CHECK(std::abs(moments.values[i] - truth) < 4.0 * moments.errors[i]);
        }
    }
    SECTION("error coverage over replicas") {
        // A standard error covers the truth in 68.3% of the replicas; with
        // 200 replicas x 4 degrees the binomial spread of the coverage is
        // 1.6%, so a window of +-5% is a 3-sigma test.
        constexpr int kReplicas = 200;
        constexpr std::size_t kSamples = 5000;
        int covered = 0;
        int trials = 0;
        for (int replica = 0; replica < kReplicas; ++replica) {
            const checks::LegendreMoments moments =
                checks::legendreMoments(sampleLegendre(kCoefficients, kSamples, rng), 4);
            for (std::size_t i = 0; i < kCoefficients.size(); ++i) {
                covered +=
                    std::abs(moments.values[i] - kCoefficients[i]) < moments.errors[i] ? 1 : 0;
                ++trials;
            }
        }
        const double coverage = static_cast<double>(covered) / trials;
        CHECK_THAT(coverage, WithinAbs(0.683, 0.05));
    }
    SECTION("isotropic sample") {
        std::uniform_real_distribution<double> cosine(-1.0, 1.0);
        std::vector<double> mu(100000);
        for (double& value : mu) {
            value = cosine(rng);
        }
        const checks::LegendreMoments moments = checks::legendreMoments(mu, kDegrees);
        for (std::size_t i = 0; i < moments.values.size(); ++i) {
            CAPTURE(i + 1);
            CHECK(std::abs(moments.values[i]) < 4.0 * moments.errors[i]);
            // Var(P_l) = 1 / (2l + 1) for a uniform cosine.
            const double expected_error =
                1.0 / std::sqrt((2.0 * static_cast<double>(i + 1) + 1.0) * 100000.0);
            CHECK_THAT(moments.errors[i], WithinRel(expected_error, 0.02));
        }
    }
    SECTION("invalid input") {
        CHECK_THROWS_AS(checks::legendreMoments({0.1, 0.2}, 0), std::invalid_argument);
        CHECK_THROWS_AS(checks::legendreMoments({0.1}, 2), std::invalid_argument);
        CHECK_THROWS_AS(checks::legendreMoments({0.1, std::nextafter(1.0, 2.0)}, 2),
                        std::invalid_argument);
        CHECK_THROWS_AS(checks::legendreMoments({0.1, kNaN}, 2), std::invalid_argument);
    }
}

namespace {

// Merging only adds small integers, which double represents exactly.
constexpr double kMergeTolerance = 1e-12;

void requireEdges(const std::vector<double>& actual, const std::vector<double>& wanted) {
    REQUIRE(actual.size() == wanted.size());
    for (std::size_t i = 0; i < wanted.size(); ++i) {
        CHECK_THAT(actual[i], WithinAbs(wanted[i], kMergeTolerance));
    }
}

} // namespace

TEST_CASE("sparse bins are merged left to right", "[checks][merge]") {
    const std::vector<double> edges{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

    SECTION("nothing to merge") {
        const auto merged = checks::mergeSparseBins({5, 6, 7, 8, 9}, {6, 6, 6, 6, 6}, edges, 5.0);
        CHECK(merged.expected.size() == 5);
        requireEdges(merged.edges, edges);
    }
    SECTION("sparse tails, totals conserved") {
        const std::vector<double> observed{1, 3, 20, 4, 0};
        const std::vector<double> expected{2, 4, 18, 3, 1};
        const auto merged = checks::mergeSparseBins(observed, expected, edges, 5.0);
        // [0,2) reaches 6; [2,3) is 18; [3,5) holds 4 < 5 and joins the last group.
        REQUIRE(merged.expected.size() == 2);
        CHECK_THAT(merged.expected[0], WithinAbs(6.0, kMergeTolerance));
        CHECK_THAT(merged.expected[1], WithinAbs(22.0, kMergeTolerance));
        CHECK_THAT(merged.observed[0], WithinAbs(4.0, kMergeTolerance));
        CHECK_THAT(merged.observed[1], WithinAbs(24.0, kMergeTolerance));
        requireEdges(merged.edges, {0.0, 2.0, 5.0});
    }
    SECTION("everything below threshold, or malformed input") {
        CHECK_THROWS_AS(
            checks::mergeSparseBins({1, 1, 1, 1, 1}, {0.5, 0.5, 0.5, 0.5, 0.5}, edges, 5.0),
            std::invalid_argument);
        CHECK_THROWS_AS(checks::mergeSparseBins({1, 1}, {6, 6}, edges, 5.0), std::invalid_argument);
        CHECK_THROWS_AS(checks::mergeSparseBins({}, {}, {0.0}, 5.0), std::invalid_argument);
        CHECK_THROWS_AS(checks::mergeSparseBins({1}, {6}, {0.0, 1.0}, 0.0), std::invalid_argument);
    }
}
