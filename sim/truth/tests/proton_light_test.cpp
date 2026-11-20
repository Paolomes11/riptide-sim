#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/sim/proton_light.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::sim::ProtonLightError;
using riptide::sim::ProtonLightTable;
namespace u = riptide::units;

constexpr double kRoundingTol = 1e-12;
constexpr double kQuadratureTol = 1e-6;
constexpr double kB = 0.155 * u::mm / u::MeV;

// Log-spaced energies from lo to hi, n of them.
std::vector<double> logGrid(double lo, double hi, std::size_t n) {
    std::vector<double> K(n);
    for (std::size_t i = 0; i < n; ++i) {
        K[i] = lo * std::pow(hi / lo, static_cast<double>(i) / static_cast<double>(n - 1));
    }
    return K;
}

// S = a / E, a power law that log-log interpolation reproduces exactly, with
// a closed-form light: integral of dE / (1 + c / E) = E - c ln(E + c), c = kB a.
constexpr double kPowerLawA = 10.0 * u::MeV * u::MeV / u::mm;
constexpr double kGridLow = 1.0 * u::keV;
constexpr double kGridHigh = 100.0 * u::MeV;

ProtonLightTable powerLawTable(double kb) {
    auto K = logGrid(kGridLow, kGridHigh, 251);
    std::vector<double> S;
    S.reserve(K.size());
    for (double E : K) {
        S.push_back(kPowerLawA / E);
    }
    return {K, S, kb};
}

double powerLawLight(double K, double kb) {
    const double c = kb * kPowerLawA;
    const double below = kGridLow / (1.0 + c / kGridLow); // S held constant below the grid
    const auto primitive = [c](double E) { return E - c * std::log(E + c); };
    return below + primitive(K) - primitive(kGridLow);
}

} // namespace

TEST_CASE("proton light: with kB = 0 the light is the energy", "[proton_light]") {
    const auto table = powerLawTable(0.0);
    CHECK_THAT(table.light(14.0 * u::MeV), WithinRel(14.0 * u::MeV, kRoundingTol));
    CHECK_THAT(table.visibleEnergy(3.0 * u::MeV, 0.5 * u::MeV),
               WithinRel(0.5 * u::MeV, kRoundingTol));
    // Deposit above the kinetic energy (nuclear-stopping overshoot): kept whole.
    CHECK_THAT(table.visibleEnergy(1.0 * u::keV, 1.2 * u::keV),
               WithinRel(1.2 * u::keV, kRoundingTol));
}

TEST_CASE("proton light: constant stopping power gives linear quenching", "[proton_light]") {
    const double S = 5.0 * u::MeV / u::mm;
    const ProtonLightTable table({1.0 * u::keV, 1.0 * u::MeV, 20.0 * u::MeV}, {S, S, S}, kB);
    const double K = 7.3 * u::MeV;
    CHECK_THAT(table.light(K), WithinRel(K / (1.0 + kB * S), kRoundingTol));
}

TEST_CASE("proton light: power-law stopping power matches the closed form", "[proton_light]") {
    const auto table = powerLawTable(kB);
    for (double K : {0.5 * u::keV, 3.0 * u::keV, 0.1 * u::MeV, 2.0 * u::MeV, 14.0 * u::MeV,
                     99.0 * u::MeV, kGridHigh}) {
        INFO("K = " << K);
        CHECK_THAT(table.light(K), WithinRel(K < kGridLow ? K / (1.0 + kB * kPowerLawA / kGridLow)
                                                          : powerLawLight(K, kB),
                                             kQuadratureTol));
    }
}

TEST_CASE("proton light: the step light telescopes to L(Ep) for any steps", "[proton_light]") {
    const auto table = powerLawTable(kB);
    const double Ep = 14.0 * u::MeV;
    // Two very different step sequences, the last step reaching zero.
    for (double fraction : {0.5, 0.01}) {
        INFO("fraction = " << fraction);
        double K = Ep;
        double light = 0.0;
        while (K > 0.0) {
            const double edep = (K < 1.0 * u::keV) ? K : fraction * K;
            light += table.visibleEnergy(K, edep);
            K -= edep;
        }
        CHECK_THAT(light, WithinRel(table.light(Ep), kRoundingTol));
    }
}

TEST_CASE("proton light: a step gives between zero and its deposit", "[proton_light]") {
    const auto table = powerLawTable(kB);
    for (double K : {0.5 * u::keV, 0.2 * u::MeV, 14.0 * u::MeV}) {
        for (double edep : {0.0, 1e-9 * u::MeV, 0.3 * K, K, 1.5 * K}) {
            INFO("K = " << K << ", edep = " << edep);
            const double evis = table.visibleEnergy(K, edep);
            CHECK(evis >= 0.0);
            CHECK(evis <= edep);
        }
    }
    CHECK_THAT(table.visibleEnergy(1.0 * u::MeV, 0.0), WithinAbs(0.0, 0.0));
}

TEST_CASE("proton light: invalid input throws", "[proton_light]") {
    using V = std::vector<double>;
    CHECK_THROWS_AS(ProtonLightTable(V{1.0}, V{1.0}, kB), ProtonLightError);
    CHECK_THROWS_AS(ProtonLightTable(V{1.0, 2.0}, V{1.0}, kB), ProtonLightError);
    CHECK_THROWS_AS(ProtonLightTable(V{2.0, 1.0}, V{1.0, 1.0}, kB), ProtonLightError);
    CHECK_THROWS_AS(ProtonLightTable(V{1.0, 2.0}, V{1.0, 0.0}, kB), ProtonLightError);
    CHECK_THROWS_AS(ProtonLightTable(V{1.0, 2.0}, V{1.0, 1.0}, -kB), ProtonLightError);
    CHECK_THROWS_AS(ProtonLightTable(V{1.0, 2.0}, V{1.0, NAN}, kB), ProtonLightError);
    const auto table = powerLawTable(kB);
    CHECK_THROWS_AS(table.light(-1.0), ProtonLightError);
    CHECK_THROWS_AS(table.light(2.0 * kGridHigh), ProtonLightError);
    CHECK_THROWS_AS(table.visibleEnergy(2.0 * kGridHigh, 1.0), ProtonLightError);
    CHECK_THROWS_AS(table.visibleEnergy(1.0, -1.0), ProtonLightError);
}
