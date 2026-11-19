#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "riptide/common/units.hpp"

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

TEST_CASE("base units are one", "[units]") {
    namespace u = riptide::units;
    REQUIRE_THAT(u::MeV, WithinULP(1.0, 0));
    REQUIRE_THAT(u::mm, WithinULP(1.0, 0));
    REQUIRE_THAT(u::ns, WithinULP(1.0, 0));
    REQUIRE_THAT(u::rad, WithinULP(1.0, 0));
}

TEST_CASE("derived units convert to base units", "[units]") {
    namespace u = riptide::units;
    REQUIRE_THAT(1500.0 * u::keV, WithinRel(1.5 * u::MeV, 1e-15));
    REQUIRE_THAT(0.02 * u::GeV, WithinRel(20.0 * u::MeV, 1e-15));
    REQUIRE_THAT(2.5 * u::cm, WithinRel(25.0 * u::mm, 1e-15));
    REQUIRE_THAT(1.2 * u::m, WithinRel(1200.0 * u::mm, 1e-15));
    REQUIRE_THAT(2.0e7 * u::eV, WithinRel(20.0 * u::MeV, 1e-15));
    REQUIRE_THAT(u::m2, WithinRel(1.0e6 * u::mm2, 1e-15));
    // 1 b = 1e-24 cm^2 = 1e-22 mm^2.
    REQUIRE_THAT(u::barn, WithinRel(1.0e-24 * u::cm * u::cm, 1e-15));
    REQUIRE_THAT(180.0 * u::deg, WithinRel(u::pi * u::rad, 1e-15));
}
