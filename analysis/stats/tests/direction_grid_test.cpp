#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/direction_grid.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using riptide::geometry::angleBetween;
using riptide::geometry::Direction;
using riptide::geometry::dot;
using riptide::stats::directionCapGrid;
using riptide::stats::directionCapRadius;
using riptide::stats::kSphereGridPixels;
using riptide::stats::sphereGrid;
using riptide::units::pi;

constexpr double kUnitTolerance = 1e-14;

// An axis with no symmetry relative to the lattice, so that no count is
// helped by an alignment.
Direction obliqueAxis() {
    return riptide::geometry::normalized(Direction{3.0, 2.0, 1.0});
}

std::size_t countWithin(const std::vector<Direction>& grid, const Direction& axis, double radius) {
    return static_cast<std::size_t>(std::count_if(
        grid.begin(), grid.end(), [&](const auto& d) { return angleBetween(d, axis) < radius; }));
}

} // namespace

TEST_CASE("sphereGrid gives unit vectors", "[direction_grid]") {
    for (const auto& d : sphereGrid(kSphereGridPixels)) {
        CHECK_THAT(dot(d, d), WithinAbs(1.0, kUnitTolerance));
    }
}

TEST_CASE("sphereGrid has equal-area pixels", "[direction_grid]") {
    // A cap of half-opening r holds the fraction (1 - cos r) / 2 of the
    // sphere. A Fibonacci lattice places points with a discrepancy far below
    // Poisson counting, so a 3-sigma Poisson band is a loose bound.
    const auto grid = sphereGrid(kSphereGridPixels);
    const auto n = static_cast<double>(grid.size());
    for (const double radius : {0.2, 0.5, 1.0, 2.0}) {
        const double expected = n * (1.0 - std::cos(radius)) / 2.0;
        const auto count = static_cast<double>(countWithin(grid, obliqueAxis(), radius));
        CHECK(std::abs(count - expected) < 3.0 * std::sqrt(expected));
    }
    // The centroid of an equal-area lattice is the origin.
    Direction sum{};
    for (const auto& d : grid) {
        sum = {sum.x + d.x, sum.y + d.y, sum.z + d.z};
    }
    CHECK(std::sqrt(dot(sum, sum)) / n < 1e-3);
}

TEST_CASE("directionCapGrid fills the cap uniformly in solid angle", "[direction_grid]") {
    const Direction centre = obliqueAxis();
    constexpr double kRadius = 0.3;
    constexpr std::size_t kPixels = 8000;
    const auto cap = directionCapGrid(centre, kRadius, kPixels);
    REQUIRE(cap.size() == kPixels);
    // cos alpha is stratified on [cos r, 1]: the empirical CDF is off the
    // uniform one by at most half a pixel, so the KS distance is <= 1/(2n)
    // up to rounding.
    std::vector<double> u;
    u.reserve(kPixels);
    for (const auto& d : cap) {
        CHECK_THAT(dot(d, d), WithinAbs(1.0, kUnitTolerance));
        const double alpha = angleBetween(d, centre);
        CHECK(alpha <= kRadius + 1e-12);
        u.push_back((1.0 - std::cos(alpha)) / (1.0 - std::cos(kRadius)));
    }
    std::sort(u.begin(), u.end());
    double ks = 0.0;
    const auto n = static_cast<double>(kPixels);
    for (std::size_t i = 0; i < u.size(); ++i) {
        const auto index = static_cast<double>(i);
        ks = std::max({ks, std::abs(u[i] - index / n), std::abs(u[i] - (index + 1.0) / n)});
    }
    CHECK(ks < 1.0 / (2.0 * n) + 1e-9);
    // The golden-angle azimuths cancel: the mean is along the centre.
    Direction sum{};
    for (const auto& d : cap) {
        sum = {sum.x + d.x, sum.y + d.y, sum.z + d.z};
    }
    CHECK(angleBetween(riptide::geometry::normalized(sum), centre) < 1e-4);
}

TEST_CASE("directionCapGrid rejects invalid input", "[direction_grid]") {
    const Direction z{0.0, 0.0, 1.0};
    CHECK_THROWS_AS(directionCapGrid(z, 0.0, 10), std::invalid_argument);
    CHECK_THROWS_AS(directionCapGrid(z, 4.0, 10), std::invalid_argument);
    CHECK_THROWS_AS(directionCapGrid(z, 0.1, 0), std::invalid_argument);
    CHECK_NOTHROW(directionCapGrid(z, pi, 10));
    CHECK_THROWS_AS(sphereGrid(0), std::invalid_argument);
}

TEST_CASE("directionCapRadius keeps the region above the threshold", "[direction_grid]") {
    const auto grid = sphereGrid(kSphereGridPixels);
    const Direction truth = obliqueAxis();
    // A Gaussian log-posterior of width 0.1 rad about the axis: pixels above
    // max - 10 lie within ~sqrt(20) * 0.1 rad of the MAP.
    constexpr double kWidth = 0.1;
    constexpr double kDeltaLog = 10.0;
    std::vector<double> log_post(grid.size());
    for (std::size_t i = 0; i < grid.size(); ++i) {
        const double a = angleBetween(grid[i], truth) / kWidth;
        log_post[i] = -0.5 * a * a;
    }
    const double radius = directionCapRadius(log_post, grid, kDeltaLog);
    const auto map_index = static_cast<std::size_t>(
        std::distance(log_post.begin(), std::max_element(log_post.begin(), log_post.end())));
    const double max_value = log_post[map_index];
    for (std::size_t i = 0; i < grid.size(); ++i) {
        if (log_post[i] > max_value - kDeltaLog) {
            CHECK(angleBetween(grid[i], grid[map_index]) < radius);
        }
    }
    // Two coarse steps on top of sqrt(2 * 10) * width; a coarse step is about
    // sqrt(4 pi / n) ~ 0.065 rad.
    CHECK(radius > std::sqrt(2.0 * kDeltaLog) * kWidth);
    CHECK(radius < std::sqrt(2.0 * kDeltaLog) * kWidth + 0.3);

    // A flat posterior keeps the whole sphere: the radius is capped at pi.
    const std::vector<double> flat(grid.size(), 0.0);
    CHECK_THAT(directionCapRadius(flat, grid, kDeltaLog), WithinAbs(pi, 0.0));
}

TEST_CASE("directionCapRadius rejects invalid input", "[direction_grid]") {
    const auto grid = sphereGrid(10);
    CHECK_THROWS_AS(directionCapRadius(std::vector<double>(9, 0.0), grid), std::invalid_argument);
    CHECK_THROWS_AS(directionCapRadius(std::vector<double>(1, 0.0), sphereGrid(1)),
                    std::invalid_argument);
    std::vector<double> bad(10, 0.0);
    bad[3] = std::nan("");
    CHECK_THROWS_AS(directionCapRadius(bad, grid), std::invalid_argument);
    CHECK_THROWS_AS(directionCapRadius(std::vector<double>(10, -INFINITY), grid),
                    std::invalid_argument);
    CHECK_THROWS_AS(directionCapRadius(std::vector<double>(10, 0.0), grid, 0.0),
                    std::invalid_argument);
}
