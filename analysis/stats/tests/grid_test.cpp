#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>
#include <stdexcept>

#include "riptide/common/units.hpp"
#include "riptide/stats/grid.hpp"

namespace st = riptide::stats;
namespace u = riptide::units;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

TEST_CASE("a grid has linspace points with exact end points", "[stats][grid]") {
    const st::Grid grid(0.5, 6.0, 500);
    REQUIRE(grid.size() == 500);
    CHECK_THAT(grid.front(), WithinULP(0.5, 0));
    CHECK_THAT(grid.back(), WithinULP(6.0, 0));
    CHECK_THAT(grid.step(), WithinRel(5.5 / 499.0, 1e-15));
    for (std::size_t i = 1; i < grid.size(); ++i) {
        CHECK_THAT(grid[i] - grid[i - 1], WithinAbs(grid.step(), 1e-14));
    }
}

TEST_CASE("invalid grids are rejected", "[stats][grid]") {
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(st::Grid(1.0, 1.0, 10), std::invalid_argument);
    CHECK_THROWS_AS(st::Grid(2.0, 1.0, 10), std::invalid_argument);
    CHECK_THROWS_AS(st::Grid(kNaN, 1.0, 10), std::invalid_argument);
    CHECK_THROWS_AS(st::Grid(0.0, std::numeric_limits<double>::infinity(), 10),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::Grid(0.0, 1.0, 1), std::invalid_argument);
}

TEST_CASE("the default grids follow D32", "[stats][grid]") {
    const st::Grid En = st::energyGrid(st::EnergyDomain{});
    CHECK(En.size() == st::kEnergyGridPoints);
    CHECK_THAT(En.front(), WithinULP(10.0, 0));
    CHECK_THAT(En.back(), WithinULP(20.0, 0));

    const st::Grid theta = st::recoilAngleGrid();
    CHECK(theta.size() == st::kRecoilAngleGridPoints);
    CHECK_THAT(theta.front(), WithinULP(0.0, 0));
    CHECK_THAT(theta.back(), WithinULP(0.5 * u::pi, 0));
}

TEST_CASE("a widened domain keeps the grid step", "[stats][grid]") {
    const st::EnergyDomain base{};
    const st::EnergyDomain wide = st::sameStepDomain(base, 30.0);
    CHECK(wide.points == 999); // 499 intervals per 10 MeV, 20 MeV wide, plus one
    CHECK_THAT(st::energyGrid(wide).step(), WithinRel(st::energyGrid(base).step(), 1e-12));
    CHECK_THAT(wide.En_min, WithinULP(base.En_min, 0));
    CHECK_THROWS_AS(st::sameStepDomain(base, 5.0), std::invalid_argument);
}
