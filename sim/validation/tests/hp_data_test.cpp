#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sstream>
#include <string>

#include "riptide/common/units.hpp"
#include "riptide/sim/hp_data.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace riptide::sim::checks;
namespace units = riptide::units;

// Values are exact decimals or one lin-lin step on them: a few ulps.
constexpr double kRel = 1.0e-14;

// Same layout as the G4NDL files, tabs and trailing whitespace included.
constexpr const char* kSample = " \t0\n \t0\n \t3\n"
                                "4.812000e+06 \t0.000000e+00 \t5.000000e+06 \t2.000000e-02 \t\n"
                                "1.400000e+07 \t4.829700e-01 \t\n";

HpCrossSection sample() {
    std::istringstream in(kSample);
    return parseHpCrossSection(in);
}

} // namespace

TEST_CASE("parseHpCrossSection reads the G4NDL layout in project units", "[hp_data]") {
    const HpCrossSection t = sample();
    REQUIRE(t.En.size() == 3);
    REQUIRE(t.sigma.size() == 3);
    CHECK_THAT(t.En[0], WithinRel(4.812 * units::MeV, kRel));
    CHECK_THAT(t.En[2], WithinRel(14.0 * units::MeV, kRel));
    CHECK_THAT(t.sigma[1], WithinRel(0.02 * units::barn, kRel));
    CHECK_THAT(t.sigma[2], WithinRel(0.48297 * units::barn, kRel));
}

TEST_CASE("evaluate interpolates lin-lin and stops at the table ends", "[hp_data]") {
    const HpCrossSection t = sample();
    CHECK_THAT(evaluate(t, 5.0 * units::MeV), WithinRel(0.02 * units::barn, kRel));
    CHECK_THAT(evaluate(t, 14.0 * units::MeV), WithinRel(0.48297 * units::barn, kRel));
    // Midpoint of the last interval.
    CHECK_THAT(evaluate(t, 9.5 * units::MeV),
               WithinRel(0.5 * (0.02 + 0.48297) * units::barn, kRel));
    CHECK_THAT(evaluate(t, t.En.front()), WithinAbs(0.0, 0.0));
    CHECK_THROWS_AS(evaluate(t, 4.0 * units::MeV), HpDataError);
    CHECK_THROWS_AS(evaluate(t, 14.5 * units::MeV), HpDataError);
}

TEST_CASE("parseHpCrossSection refuses malformed files", "[hp_data]") {
    for (const std::string text :
         {"", "0 0", "0 0 0", "0 0 -2", "0 0 2 1e6 1.0", "0 0 3 1e6 1.0 1e6 2.0 1e6 3.0",
          "0 0 2 2e6 1.0 1e6 2.0", "0 0 1 1e6 -1.0", "0 0 1 1e6 abc"}) {
        CAPTURE(text);
        std::istringstream in(text);
        CHECK_THROWS_AS(parseHpCrossSection(in), HpDataError);
    }
}

TEST_CASE("a repeated energy is a discontinuity, valued above it", "[hp_data]") {
    std::istringstream in("0 0 4 1e6 1.0 2e6 2.0 2e6 5.0 3e6 6.0");
    const HpCrossSection t = parseHpCrossSection(in);
    REQUIRE(t.En.size() == 4);
    CHECK_THAT(evaluate(t, 2.0 * units::MeV), WithinRel(5.0 * units::barn, kRel));
    CHECK_THAT(evaluate(t, 1.5 * units::MeV), WithinRel(1.5 * units::barn, kRel));
    CHECK_THAT(evaluate(t, 2.5 * units::MeV), WithinRel(5.5 * units::barn, kRel));
}
