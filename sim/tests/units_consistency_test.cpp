// Values cross from Geant4 to the analysis without conversion only if the
// project units coincide exactly with Geant4's internal (CLHEP) units.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <CLHEP/Units/SystemOfUnits.h>

#include "riptide/common/units.hpp"

using Catch::Matchers::WithinULP;

TEST_CASE("project units equal Geant4 internal units", "[units][geant4]") {
    namespace u = riptide::units;
    CHECK_THAT(u::MeV, WithinULP(CLHEP::MeV, 0));
    CHECK_THAT(u::keV, WithinULP(CLHEP::keV, 0));
    CHECK_THAT(u::GeV, WithinULP(CLHEP::GeV, 0));
    CHECK_THAT(u::mm, WithinULP(CLHEP::mm, 0));
    CHECK_THAT(u::cm, WithinULP(CLHEP::cm, 0));
    CHECK_THAT(u::m, WithinULP(CLHEP::m, 0));
    CHECK_THAT(u::ns, WithinULP(CLHEP::ns, 0));
    CHECK_THAT(u::rad, WithinULP(CLHEP::rad, 0));
    CHECK_THAT(u::deg, WithinULP(CLHEP::deg, 1));
}
