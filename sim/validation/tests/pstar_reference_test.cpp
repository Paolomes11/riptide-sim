// Reads the PSTAR table on the data disk and compares a few rows with the
// numbers copied from the file by hand. The disk may be absent (CI, unplugged
// SSD): then the case is skipped.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/pstar.hpp"

using Catch::Matchers::WithinRel;

namespace checks = riptide::sim::checks;
namespace units = riptide::units;

namespace {

// Tabulated values are exact decimals times the density: a few ulps.
constexpr double kGridTolerance = 1e-13;
constexpr std::size_t kRows = 132;
constexpr double kDensity = checks::kVinyltolueneDensityGramPerCm3;
constexpr double kStopping = kDensity * units::MeV / units::cm;
constexpr double kRange = units::cm / kDensity;
constexpr double kB = 0.155; // mm/MeV, the D18 value
constexpr double kConvergence = 1e-9;

checks::PstarTable readOrSkip() {
    const std::filesystem::path path =
        riptide::data::resolveDataRoot(std::nullopt) / "reference" / "pstar" / "pstar_216.txt";
    if (!std::filesystem::is_regular_file(path)) {
        SKIP("reference file not available: pstar/pstar_216.txt");
    }
    std::ifstream in(path);
    REQUIRE(in);
    return checks::parsePstarTable(in, kDensity);
}

} // namespace

TEST_CASE("PSTAR vinyltoluene-based plastic scintillator", "[reference][pstar]") {
    const checks::PstarTable t = readOrSkip();
    REQUIRE(t.Ep.size() == kRows);
    CHECK_THAT(t.Ep.front(), WithinRel(1.0e-3 * units::MeV, kGridTolerance));
    CHECK_THAT(t.Ep.back(), WithinRel(1.0e4 * units::MeV, kGridTolerance));

    // 1.000E+01 4.577E+01 2.463E-02 4.579E+01 1.222E-01 1.220E-01 0.9984
    CHECK_THAT(checks::stoppingPower(t, 10.0), WithinRel(45.79 * kStopping, kGridTolerance));
    CHECK_THAT(checks::csdaRange(t, 10.0), WithinRel(0.1222 * kRange, kGridTolerance));
    // 2.000E+00 1.605E+02 1.102E-01 1.606E+02 7.388E-03 7.363E-03 0.9967
    CHECK_THAT(checks::stoppingPower(t, 2.0), WithinRel(160.6 * kStopping, kGridTolerance));
    // 2.000E+01 2.603E+01 1.280E-02 2.604E+01 4.250E-01 4.244E-01 0.9987
    CHECK_THAT(checks::csdaRange(t, 20.0), WithinRel(0.4250 * kRange, kGridTolerance));

    // 14 MeV is between the 12.5 and 15 MeV rows: the interpolated range must
    // lie between theirs, and the light below the energy.
    CHECK(checks::csdaRange(t, 14.0) > checks::csdaRange(t, 12.5));
    CHECK(checks::csdaRange(t, 14.0) < checks::csdaRange(t, 15.0));
    const double light = checks::visibleEnergy(t, 14.0, kB);
    CHECK(light > 0.0);
    CHECK(light < 14.0);
}

TEST_CASE("PSTAR light: the default Simpson panels have converged on the real table",
          "[reference][pstar]") {
    const checks::PstarTable t = readOrSkip();
    for (const double Ep : {0.37, 2.0, 10.0, 14.0, 20.0}) {
        CAPTURE(Ep);
        CHECK_THAT(
            checks::visibleEnergy(t, Ep, kB),
            WithinRel(checks::visibleEnergy(t, Ep, kB, 2 * checks::kBirksPanels), kConvergence));
    }
}
