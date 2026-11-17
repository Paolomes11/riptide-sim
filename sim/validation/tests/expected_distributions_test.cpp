#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/expected_distributions.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace checks = riptide::sim::checks;
namespace units = riptide::units;
using riptide::schema::FirstInteraction;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kPerCubicMm = 1.0 / (units::mm * units::mm * units::mm);

// Coefficients of a markedly anisotropic, positive distribution.
const std::vector<double> kCoefficients = {0.3, 0.1, -0.05, 0.02};

double legendreDensity(const std::vector<double>& a, double mu) {
    double value = 0.5;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const int degree = static_cast<int>(i) + 1;
        value += 0.5 * (2.0 * degree + 1.0) * a[i] * riptide::kinematics::legendreP(degree, mu);
    }
    return value;
}

// Composite Simpson rule of f over [-1, mu]: an independent check of the
// closed-form integral, exact here up to O(h^4) on a degree-4 polynomial.
double simpsonCdf(const std::vector<double>& a, double mu) {
    constexpr int kIntervals = 2000; // even
    const double h = (mu + 1.0) / kIntervals;
    double sum = legendreDensity(a, -1.0) + legendreDensity(a, mu);
    for (int i = 1; i < kIntervals; ++i) {
        sum += (i % 2 == 1 ? 4.0 : 2.0) * legendreDensity(a, -1.0 + i * h);
    }
    return sum * h / 3.0;
}

double sum(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0);
}

std::vector<double> uniformEdges(double low, double high, int bins) {
    std::vector<double> edges;
    for (int i = 0; i <= bins; ++i) {
        edges.push_back(low + (high - low) * i / bins);
    }
    return edges;
}

} // namespace

TEST_CASE("Legendre CDF", "[expected][legendre]") {
    CHECK_THAT(checks::legendreCdf(kCoefficients, -1.0), WithinAbs(0.0, 1e-15));
    CHECK_THAT(checks::legendreCdf(kCoefficients, 1.0), WithinRel(1.0, 1e-15));
    for (const double mu : {-0.9, -0.5, 0.0, 0.3, 0.77, 0.99}) {
        CAPTURE(mu);
        CHECK_THAT(checks::legendreCdf(kCoefficients, mu),
                   WithinAbs(simpsonCdf(kCoefficients, mu), 1e-12));
        // One coefficient: F = (mu + 1)/2 + 3/4 a_1 (mu^2 - 1).
        CHECK_THAT(checks::legendreCdf({0.2}, mu),
                   WithinAbs(0.5 * (mu + 1.0) + 0.75 * 0.2 * (mu * mu - 1.0), 1e-15));
        // Isotropy.
        CHECK_THAT(checks::legendreCdf({}, mu), WithinAbs(0.5 * (mu + 1.0), 1e-15));
    }
    CHECK_THROWS_AS(checks::legendreCdf(kCoefficients, 1.0 + 1e-12), std::invalid_argument);
    CHECK_THROWS_AS(checks::legendreCdf(kCoefficients, kNaN), std::invalid_argument);
}

TEST_CASE("proton energy spectrum from Legendre coefficients", "[expected][legendre]") {
    const std::vector<double> edges = uniformEdges(0.0, 1.0, 20);
    SECTION("isotropy gives a flat spectrum") {
        for (const double p : checks::protonEnergyBinProbabilities({}, edges)) {
            CHECK_THAT(p, WithinRel(0.05, 1e-13));
        }
    }
    SECTION("anisotropic: normalized, and each bin is the cosine interval") {
        const std::vector<double> p = checks::protonEnergyBinProbabilities(kCoefficients, edges);
        REQUIRE(p.size() == 20);
        CHECK_THAT(sum(p), WithinRel(1.0, 1e-14));
        // Low proton energy = forward neutron: the first bin is mu in [0.9, 1].
        CHECK_THAT(
            p.front(),
            WithinAbs(simpsonCdf(kCoefficients, 1.0) - simpsonCdf(kCoefficients, 0.9), 1e-12));
        // a_1 > 0 favours forward neutrons, hence low proton energies.
        CHECK(p.front() > p.back());
    }
    SECTION("a partial range does not sum to 1") {
        CHECK(sum(checks::protonEnergyBinProbabilities(kCoefficients, {0.2, 0.5, 0.8})) < 1.0);
    }
    SECTION("invalid edges") {
        CHECK_THROWS_AS(checks::protonEnergyBinProbabilities({}, {0.5}), std::invalid_argument);
        CHECK_THROWS_AS(checks::protonEnergyBinProbabilities({}, {0.0, 1.1}),
                        std::invalid_argument);
        CHECK_THROWS_AS(checks::protonEnergyBinProbabilities({}, {0.0, 0.5, 0.5}),
                        std::invalid_argument);
    }
}

TEST_CASE("depth of the first interaction", "[expected][depth]") {
    constexpr double kThickness = 50.0 * units::mm;
    const std::vector<double> edges = uniformEdges(0.0, kThickness, 10);
    SECTION("normalized and decreasing") {
        const std::vector<double> p =
            checks::truncatedExponentialBinProbabilities(0.02 / units::mm, kThickness, edges);
        CHECK_THAT(sum(p), WithinRel(1.0, 1e-14));
        // Ratio of consecutive bins is e^{-Sigma dz}.
        for (std::size_t i = 1; i < p.size(); ++i) {
            CHECK_THAT(p[i] / p[i - 1], WithinRel(std::exp(-0.02 * 5.0), 1e-12));
        }
    }
    SECTION("vanishing Sigma tends to the uniform distribution") {
        for (const double p :
             checks::truncatedExponentialBinProbabilities(0.0, kThickness, edges)) {
            CHECK_THAT(p, WithinRel(0.1, 1e-14));
        }
        // A tiny Sigma must not lose its digits to 1 - e^{-Sigma L}.
        for (const double p :
             checks::truncatedExponentialBinProbabilities(1e-12 / units::mm, kThickness, edges)) {
            CHECK_THAT(p, WithinRel(0.1, 1e-9));
        }
    }
    SECTION("invalid input") {
        CHECK_THROWS_AS(checks::truncatedExponentialBinProbabilities(-1.0, kThickness, edges),
                        std::invalid_argument);
        CHECK_THROWS_AS(checks::truncatedExponentialBinProbabilities(0.02, 0.0, edges),
                        std::invalid_argument);
        CHECK_THROWS_AS(checks::truncatedExponentialBinProbabilities(0.02, 40.0, edges),
                        std::invalid_argument);
    }
}

TEST_CASE("first-interaction probabilities from cross sections", "[expected][first]") {
    // Round numbers, not physical values: the test checks the bookkeeping.
    const checks::MediumComponent hydrogen{checks::Target::Hydrogen, 5.0e19 * kPerCubicMm,
                                           0.9 * units::barn, 0.0001 * units::barn,
                                           0.9001 * units::barn};
    const checks::MediumComponent carbon{checks::Target::Carbon, 4.0e19 * kPerCubicMm,
                                         0.7 * units::barn, 0.0002 * units::barn,
                                         1.2 * units::barn};
    constexpr double kLength = 50.0 * units::mm;

    SECTION("hydrogen and carbon") {
        const checks::FirstInteractionProbabilities p =
            checks::firstInteraction({hydrogen, carbon}, kLength);
        const double sigma_h = 5.0e19 * 0.9001e-22;
        const double sigma_c = 4.0e19 * 1.2e-22;
        const double total = sigma_h + sigma_c; // 1/mm, with 1 b = 1e-22 mm^2
        CHECK_THAT(p.sigma_total, WithinRel(total, 1e-12));
        const double interacts = 1.0 - std::exp(-total * kLength);
        CHECK_THAT(p.probability(FirstInteraction::Unknown),
                   WithinRel(std::exp(-total * kLength), 1e-12));
        CHECK_THAT(p.probabilityGivenInteraction(FirstInteraction::NpElastic),
                   WithinRel(5.0e19 * 0.9e-22 / total, 1e-12));
        CHECK_THAT(p.probabilityGivenInteraction(FirstInteraction::NCElastic),
                   WithinRel(4.0e19 * 0.7e-22 / total, 1e-12));
        CHECK_THAT(p.probabilityGivenInteraction(FirstInteraction::NCInelastic),
                   WithinRel(4.0e19 * (1.2 - 0.7 - 0.0002) * 1e-22 / total, 1e-9));
        CHECK_THAT(p.probabilityGivenInteraction(FirstInteraction::Other),
                   WithinRel((5.0e19 * 0.0001 + 4.0e19 * 0.0002) * 1e-22 / total, 1e-9));
        CHECK_THAT(p.probability(FirstInteraction::NpElastic),
                   WithinRel(interacts * 5.0e19 * 0.9e-22 / total, 1e-12));
        CHECK_THAT(std::accumulate(p.unconditional.begin(), p.unconditional.end(), 0.0),
                   WithinRel(1.0, 1e-14));
        CHECK_THAT(std::accumulate(p.conditional.begin(), p.conditional.end(), 0.0),
                   WithinRel(1.0, 1e-14));
        CHECK_THAT(p.probabilityGivenInteraction(FirstInteraction::Unknown), WithinAbs(0.0, 0.0));
    }
    SECTION("one component gives sigma_el / sigma_tot") {
        const checks::FirstInteractionProbabilities p = checks::firstInteraction({carbon}, kLength);
        CHECK_THAT(p.probabilityGivenInteraction(FirstInteraction::NCElastic),
                   WithinRel(0.7 / 1.2, 1e-12));
        CHECK_THAT(p.probabilityGivenInteraction(FirstInteraction::NpElastic), WithinAbs(0.0, 0.0));
    }
    SECTION("a thin slab almost never interacts, with unchanged branching") {
        const checks::FirstInteractionProbabilities thin =
            checks::firstInteraction({hydrogen, carbon}, 1e-9 * units::mm);
        const checks::FirstInteractionProbabilities thick =
            checks::firstInteraction({hydrogen, carbon}, kLength);
        CHECK_THAT(thin.probability(FirstInteraction::Unknown), WithinRel(1.0, 1e-10));
        CHECK_THAT(thin.probability(FirstInteraction::NpElastic),
                   WithinRel(5.0e19 * 0.9e-22 * 1e-9, 1e-6));
        CHECK_THAT(
            thin.probabilityGivenInteraction(FirstInteraction::NpElastic),
            WithinRel(thick.probabilityGivenInteraction(FirstInteraction::NpElastic), 1e-14));
    }
    SECTION("invalid input") {
        checks::MediumComponent inconsistent = carbon;
        inconsistent.sigma_elastic = 1.3 * units::barn;
        CHECK_THROWS_AS(checks::firstInteraction({inconsistent}, kLength), std::invalid_argument);
        CHECK_THROWS_AS(checks::firstInteraction({}, kLength), std::invalid_argument);
        CHECK_THROWS_AS(checks::firstInteraction({carbon}, 0.0), std::invalid_argument);
        checks::MediumComponent empty = carbon;
        empty.number_density = 0.0;
        CHECK_THROWS_AS(checks::firstInteraction({empty}, kLength), std::invalid_argument);
    }
}
