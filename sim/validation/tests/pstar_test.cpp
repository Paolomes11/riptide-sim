#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/sim/pstar.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace riptide::sim::checks;
namespace units = riptide::units;

// Values are exact decimals times the density, or one log-log step on them.
constexpr double kRel = 1.0e-13;
constexpr double kDensity = kVinyltolueneDensityGramPerCm3;
constexpr double kB = 0.155; // mm/MeV, the D18 value

// Header and three rows copied from the NIST file, trailing spaces included.
constexpr const char* kSample =
    "PSTAR: Stopping Powers and Range Tables for Protons\n"
    "\n"
    "PLASTIC SCINTILLATOR (vinyltoluene-based)                               \n"
    "\n"
    "Kinetic   Electron. Nuclear   Total     CSDA      Projected Detour \n"
    "Energy    Stp. Pow. Stp. Pow. Stp. Pow. Range     Range     Factor \n"
    "MeV       MeV cm2/g MeV cm2/g MeV cm2/g g/cm2     g/cm2            \n"
    "\n"
    "1.000E+01 4.577E+01 2.463E-02 4.579E+01 1.222E-01 1.220E-01 0.9984 \n"
    "1.500E+01 3.292E+01 1.680E-02 3.294E+01 2.529E-01 2.526E-01 0.9986 \n"
    "2.000E+01 2.603E+01 1.280E-02 2.604E+01 4.250E-01 4.244E-01 0.9987 \n";

PstarTable sample() {
    std::istringstream in(kSample);
    return parsePstarTable(in, kDensity);
}

// S = amplitude * E^exponent on the given energies, exact for log-log
// interpolation; ranges set to the same law, which is enough for the tests.
PstarTable powerLaw(const std::vector<double>& energies, double amplitude, double exponent) {
    PstarTable t;
    for (const double E : energies) {
        const double S = amplitude * std::pow(E, exponent);
        t.Ep.push_back(E);
        t.S_el.push_back(S);
        t.S_nuc.push_back(0.0);
        t.S_tot.push_back(S);
        t.csda_range.push_back(S);
        t.projected_range.push_back(S);
    }
    return t;
}

// Shape of the real table above the Bragg peak: S falling as a power of E,
// of the order of the PVT values (about 5 MeV/mm at 10 MeV), on energies a
// factor 2 apart, the widest spacing of the NIST table.
PstarTable fallingTable() {
    constexpr double kFirst = 0.001;
    constexpr double kLast = 20.0;
    std::vector<double> energies;
    for (double E = kFirst; E < kLast; E *= 2.0) {
        energies.push_back(E);
    }
    energies.push_back(kLast);
    return powerLaw(energies, 30.0, -0.8);
}

} // namespace

TEST_CASE("parsePstarTable reads the NIST layout in project units", "[pstar]") {
    const PstarTable t = sample();
    REQUIRE(t.Ep.size() == 3);
    CHECK_THAT(t.Ep[0], WithinRel(10.0 * units::MeV, kRel));
    CHECK_THAT(t.Ep[2], WithinRel(20.0 * units::MeV, kRel));
    // (MeV cm2/g) * (g/cm3) = MeV/cm; (g/cm2) / (g/cm3) = cm.
    CHECK_THAT(t.S_tot[0], WithinRel(45.79 * kDensity * units::MeV / units::cm, kRel));
    CHECK_THAT(t.S_el[1], WithinRel(32.92 * kDensity * units::MeV / units::cm, kRel));
    CHECK_THAT(t.S_nuc[2], WithinRel(0.0128 * kDensity * units::MeV / units::cm, kRel));
    CHECK_THAT(t.csda_range[1], WithinRel(0.2529 / kDensity * units::cm, kRel));
    CHECK_THAT(t.projected_range[2], WithinRel(0.4244 / kDensity * units::cm, kRel));
}

TEST_CASE("parsePstarTable refuses malformed tables", "[pstar]") {
    const std::string header = "PSTAR: Stopping Powers\nMeV MeV cm2/g\n";
    const std::string good = "1.0 2.0 1.0 3.0 1.0 1.0 1.0\n";
    for (const std::string& text : {
             std::string{},                                     // empty
             header,                                            // header only
             good,                                              // no PSTAR header
             header + good + "1.0 2.0 1.0 3.0 1.0 1.0\n",       // six columns
             header + good + "2.0 2.0 1.0 3.0 1.0 1.0 1.0 9\n", // eight columns
             header + good + "1.0 2.0 1.0 3.0 1.0 1.0 1.0\n",   // repeated energy
             header + "2.0 2.0 1.0 3.0 1.0 1.0 1.0\n" + good,   // decreasing energy
             header + "1.0 2.0 1.0 3.1 1.0 1.0 1.0\n",          // total != sum
             header + "1.0 -2.0 1.0 -1.0 1.0 1.0 1.0\n",        // negative
             header + "1.0 2.0 0.0 2.0 1.0 1.0 1.0\n",          // zero
             header + good + "trailing text\n",                 // text after the data
         }) {
        CAPTURE(text);
        std::istringstream in(text);
        CHECK_THROWS_AS(parsePstarTable(in, kDensity), PstarError);
    }
    std::istringstream in(header + good);
    CHECK_THROWS_AS(parsePstarTable(in, 0.0), PstarError);
}

TEST_CASE("parsePstarTable accepts blank lines inside and after the data", "[pstar]") {
    std::istringstream in("PSTAR\n1.0 2.0 1.0 3.0 1.0 1.0 1.0\n\n2.0 2.0 1.0 3.0 1.0 1.0 1.0\n \n");
    CHECK(parsePstarTable(in, kDensity).Ep.size() == 2);
}

TEST_CASE("log-log interpolation is exact on a power law", "[pstar]") {
    const PstarTable t = powerLaw({1.0, 2.0, 4.0}, 7.0, -0.8);
    for (const double E : {1.0, 1.5, 2.0, 3.0, 4.0}) {
        CAPTURE(E);
        CHECK_THAT(stoppingPower(t, E), WithinRel(7.0 * std::pow(E, -0.8), kRel));
        CHECK_THAT(csdaRange(t, E), WithinRel(7.0 * std::pow(E, -0.8), kRel));
    }
    CHECK_THROWS_AS(stoppingPower(t, 0.5), PstarError);
    CHECK_THROWS_AS(csdaRange(t, 4.5), PstarError);
}

TEST_CASE("with kB = 0 the light is the proton energy", "[pstar]") {
    const PstarTable t = fallingTable();
    for (const double Ep : {0.0, 0.0005, 0.001, 0.37, 2.0, 14.0, 20.0}) {
        CAPTURE(Ep);
        CHECK_THAT(visibleEnergy(t, Ep, 0.0), WithinRel(Ep, kRel));
    }
}

TEST_CASE("with constant S the light has its closed form", "[pstar]") {
    const double S = 5.0 * units::MeV / units::mm;
    const PstarTable t = powerLaw({0.001, 1.0, 10.0, 20.0}, S, 0.0);
    for (const double Ep : {0.0005, 0.5, 10.0, 14.0, 20.0}) {
        CAPTURE(Ep);
        CHECK_THAT(visibleEnergy(t, Ep, kB), WithinRel(Ep / (1.0 + kB * S), kRel));
    }
}

TEST_CASE("the light grows with Ep, falls with kB and stays below Ep", "[pstar]") {
    const PstarTable t = fallingTable();
    double previous = 0.0;
    for (const double Ep : {0.01, 0.5, 2.0, 5.0, 10.0, 14.0, 20.0}) {
        CAPTURE(Ep);
        const double light = visibleEnergy(t, Ep, kB);
        CHECK(light > previous);
        CHECK(light < Ep);
        CHECK(visibleEnergy(t, Ep, kB + 0.005) < light);
        previous = light;
    }
}

TEST_CASE("Simpson's rule converges as h^4 and the default panels suffice", "[pstar]") {
    // Spacing a factor 2 everywhere is harsher than the NIST table, where the
    // reference test asks for 1e-9; here 1e-8, and the order of the rule.
    constexpr double kTolerance = 1.0e-8;
    constexpr double kMinRatio = 12.0; // 2^4 = 16 for an h^4 error
    constexpr double kMaxRatio = 20.0;
    const PstarTable t = fallingTable();
    for (const double Ep : {0.37, 2.0, 14.0, 20.0}) {
        CAPTURE(Ep);
        const double l1 = visibleEnergy(t, Ep, kB);
        const double l2 = visibleEnergy(t, Ep, kB, 2 * kBirksPanels);
        const double l4 = visibleEnergy(t, Ep, kB, 4 * kBirksPanels);
        CHECK_THAT(l1, WithinRel(l2, kTolerance));
        const double ratio = (l1 - l2) / (l2 - l4);
        CHECK(ratio > kMinRatio);
        CHECK(ratio < kMaxRatio);
    }
}

TEST_CASE("visibleEnergy refuses invalid arguments", "[pstar]") {
    const PstarTable t = fallingTable();
    CHECK_THROWS_AS(visibleEnergy(t, -1.0, kB), PstarError);
    CHECK_THROWS_AS(visibleEnergy(t, 21.0, kB), PstarError);
    CHECK_THROWS_AS(visibleEnergy(t, 14.0, -0.1), PstarError);
    CHECK_THROWS_AS(visibleEnergy(t, 14.0, kB, 0), PstarError);
    CHECK_THROWS_AS(visibleEnergy(PstarTable{}, 1.0, kB), PstarError);
    CHECK_THAT(visibleEnergy(t, 0.0, kB), WithinAbs(0.0, 0.0));
}

TEST_CASE("BirksLightTable gives the light of visibleEnergy", "[pstar]") {
    // The same sums in the same order: equal to rounding.
    constexpr double kSameSum = 1.0e-12;
    const PstarTable t = fallingTable();
    const BirksLightTable table(t, kB);
    CHECK_THAT(table.kB(), WithinAbs(kB, 0.0));
    // Below, on and between the tabulated energies, and at both ends.
    for (const double Ep : {0.0, 0.0005, 0.001, 0.002, 0.003, 0.37, 2.048, 9.9, 14.0, 20.0}) {
        CAPTURE(Ep);
        CHECK_THAT(table.light(Ep), WithinRel(visibleEnergy(t, Ep, kB), kSameSum));
    }
    const BirksLightTable unquenched(t, 0.0);
    CHECK_THAT(unquenched.light(14.0), WithinRel(14.0, kRel));
}

TEST_CASE("BirksLightTable refuses invalid arguments", "[pstar]") {
    const PstarTable t = fallingTable();
    CHECK_THROWS_AS(BirksLightTable(PstarTable{}, kB), PstarError);
    CHECK_THROWS_AS(BirksLightTable(t, -0.1), PstarError);
    CHECK_THROWS_AS(BirksLightTable(t, kB, 0), PstarError);
    const BirksLightTable table(t, kB);
    CHECK_THROWS_AS(table.light(-1.0), PstarError);
    CHECK_THROWS_AS(table.light(21.0), PstarError);
}
