#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/sim/pstar_checks.hpp"

namespace checks = riptide::sim::checks;
namespace units = riptide::units;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::sim::DepositRecord;
using riptide::sim::TruthRecord;

namespace {

constexpr double kEnMax = 20.0 * units::MeV;
constexpr double kHalfEdge = 50.0 * units::mm;
constexpr double kKB = 0.155 * units::mm / units::MeV;
constexpr std::size_t kBins = 20; // 1 MeV bins over [0, 20] MeV
constexpr double kTolerance = 1e-12;

checks::DepositSetup setup(double kB = kKB) {
    return {{kHalfEdge, kHalfEdge, kHalfEdge}, kEnMax, kB, kBins};
}

// Falling stopping power and growing range, power laws of the order of the
// plastic values, on energies a factor 2 apart up to En_max.
checks::PstarTable table() {
    checks::PstarTable t;
    std::vector<double> energies;
    for (double E = 0.001; E < kEnMax; E *= 2.0) {
        energies.push_back(E);
    }
    energies.push_back(kEnMax);
    for (const double E : energies) {
        const double S = 30.0 * std::pow(E, -0.8);
        t.Ep.push_back(E);
        t.S_el.push_back(S);
        t.S_nuc.push_back(0.0);
        t.S_tot.push_back(S);
        t.csda_range.push_back(0.02 * std::pow(E, 1.8));
        t.projected_range.push_back(0.02 * std::pow(E, 1.8));
    }
    return t;
}

// A recoil proton stopping inside the block with the given path length and light.
struct Proton {
    TruthRecord truth;
    DepositRecord deposit;
};

Proton proton(double Ep, double length, double Evis) {
    Proton p;
    p.truth.Ep = Ep;
    p.truth.valid_np = true;
    p.deposit.proton.seen = true;
    p.deposit.proton.length = length;
    p.deposit.proton.Edep = Ep;
    p.deposit.proton.Evis = Evis;
    return p;
}

// A proton whose range and light are exactly the table's, shifted by the offsets.
Proton onTable(double Ep, double range_offset = 0.0, double light_offset = 0.0, double kB = kKB) {
    const checks::PstarTable t = table();
    return proton(Ep, checks::csdaRange(t, Ep) + range_offset,
                  checks::visibleEnergy(t, Ep, kB) + light_offset);
}

checks::PstarResults compare(const std::vector<Proton>& protons, double kB = kKB) {
    checks::PstarComparison c(setup(kB), table());
    for (const Proton& p : protons) {
        c.add(p.truth, p.deposit);
    }
    return c.finish();
}

} // namespace

TEST_CASE("protons on the table leave zero residuals", "[pstar_checks]") {
    const checks::PstarResults r = compare({onTable(3.2), onTable(13.7), onTable(19.9)});
    CHECK(r.samples == 3);
    CHECK(r.below_table == 0);
    CHECK(r.range.samples() == 3);
    CHECK_THAT(r.range.residualMean(), WithinAbs(0.0, kTolerance));
    CHECK_THAT(r.light.residualMean(), WithinAbs(0.0, kTolerance));
    CHECK_THAT(*r.range.residualSd(), WithinAbs(0.0, kTolerance));
    CHECK_THAT(*r.light.relativeResidual(), WithinAbs(0.0, kTolerance));
}

TEST_CASE("only contained protons without hadronic interactions or nuclear recoils are compared",
          "[pstar_checks]") {
    Proton exits = onTable(10.0);
    exits.deposit.proton.exits = true;
    Proton hadronic = onTable(10.0);
    hadronic.deposit.proton.n_hadronic = 1;
    Proton recoil = onTable(10.0);
    recoil.deposit.proton.n_nuclear_recoil = 1;
    Proton invalid = onTable(10.0);
    invalid.truth.valid_np = false;
    const Proton below = proton(0.0005, 1e-6, 0.0004);
    const checks::PstarResults r =
        compare({exits, hadronic, recoil, invalid, below, onTable(10.0)});
    CHECK(r.samples == 1);
    CHECK(r.below_table == 1);
    CHECK(r.range.samples() == 1);
    CHECK(r.light.samples() == 1);
}

TEST_CASE("each proton goes to the bin of its Ep, edges as in the deposit checks",
          "[pstar_checks]") {
    const checks::PstarResults r = compare({onTable(13.7), onTable(13.2), onTable(2.5)});
    REQUIRE(r.edges.size() == kBins + 1);
    CHECK_THAT(r.edges.front(), WithinAbs(0.0, 0.0));
    CHECK_THAT(r.edges.back(), WithinAbs(kEnMax, 0.0));
    REQUIRE(r.range_bins.size() == kBins);
    REQUIRE(r.light_bins.size() == kBins);
    CHECK(r.range_bins.at(13).samples() == 2);
    CHECK(r.light_bins.at(2).samples() == 1);
    CHECK(r.range_bins.at(12).samples() == 0);
    // The expectation is the mean over the protons, not the bin centre.
    const checks::PstarTable t = table();
    const double expected = 0.5 * (checks::csdaRange(t, 13.7) + checks::csdaRange(t, 13.2));
    CHECK_THAT(r.range_bins.at(13).expectedMean(), WithinRel(expected, kTolerance));
    // A single proton has a mean but no spread.
    CHECK_FALSE(r.light_bins.at(2).residualSd().has_value());
    CHECK_FALSE(r.light_bins.at(2).relativeSem().has_value());
    CHECK_FALSE(r.light_bins.at(12).relativeResidual().has_value());
}

TEST_CASE("a known offset and two values give the closed-form statistics", "[pstar_checks]") {
    constexpr double kRangeOffset = 0.03 * units::mm;
    constexpr double kLightOffset = -0.02 * units::MeV;
    const checks::PstarResults shifted = compare(
        {onTable(11.0, kRangeOffset, kLightOffset), onTable(17.0, kRangeOffset, kLightOffset)});
    CHECK_THAT(shifted.range.residualMean(), WithinRel(kRangeOffset, kTolerance));
    CHECK_THAT(shifted.light.residualMean(), WithinRel(kLightOffset, kTolerance));
    CHECK_THAT(*shifted.range.residualSd(), WithinAbs(0.0, kTolerance));
    CHECK_THAT(*shifted.light.relativeResidual(),
               WithinRel(kLightOffset / shifted.light.expectedMean(), kTolerance));

    // d = +a and -a: mean 0, sample sd a * sqrt(2), standard error a.
    constexpr double a = 0.05;
    const checks::PstarResults pair = compare({onTable(12.0, a, a), onTable(12.0, -a, -a)});
    const checks::ResidualStats& range = pair.range_bins.at(12);
    CHECK_THAT(range.residualMean(), WithinAbs(0.0, kTolerance));
    CHECK_THAT(*range.residualSd(), WithinRel(a * std::sqrt(2.0), kTolerance));
    CHECK_THAT(*range.residualSem(), WithinRel(a, kTolerance));
    CHECK_THAT(*range.relativeSem(), WithinRel(a / range.expectedMean(), kTolerance));
    CHECK_THAT(range.observedMean(), WithinRel(range.expectedMean(), kTolerance));
}

TEST_CASE("with kB = 0 the expected light is the proton energy", "[pstar_checks]") {
    const checks::PstarResults r = compare({proton(14.0, 2.0, 14.0)}, 0.0);
    CHECK_THAT(r.light.expectedMean(), WithinRel(14.0, kTolerance));
    CHECK_THAT(r.light.residualMean(), WithinAbs(0.0, kTolerance));
}

TEST_CASE("ResidualStats keeps the spread of values far from zero", "[pstar_checks]") {
    // A naive sum of squares would lose the 0.5 spread next to 1e8^2.
    constexpr double kLarge = 1.0e8;
    checks::ResidualStats s;
    for (int i = 0; i < 1000; ++i) {
        s.add(kLarge + ((i % 2 == 0) ? 0.5 : -0.5), 0.0);
    }
    CHECK_THAT(s.residualMean(), WithinRel(kLarge, kTolerance));
    // Sample sd of n values +-0.5: 0.5 * sqrt(n / (n - 1)).
    CHECK_THAT(*s.residualSd(), WithinRel(0.5 * std::sqrt(1000.0 / 999.0), 1e-9));
    CHECK_FALSE(checks::ResidualStats{}.relativeResidual().has_value());
}

TEST_CASE("PstarComparison refuses an invalid setup", "[pstar_checks]") {
    checks::DepositSetup s = setup();
    s.bins = 0;
    CHECK_THROWS_AS(checks::PstarComparison(s, table()), std::invalid_argument);
    s = setup();
    s.En_max = 0.0;
    CHECK_THROWS_AS(checks::PstarComparison(s, table()), std::invalid_argument);
    s = setup();
    s.En_max = 25.0; // above the table
    CHECK_THROWS_AS(checks::PstarComparison(s, table()), checks::PstarError);
    CHECK_THROWS_AS(checks::PstarComparison(setup(-0.1), table()), checks::PstarError);
}
