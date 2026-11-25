#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "riptide/common/direction.hpp"
#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using riptide::geometry::angleBetween;
using riptide::geometry::cross;
using riptide::geometry::Direction;
using riptide::geometry::directionAroundAxis;
using riptide::geometry::dot;
using riptide::geometry::perpendicularBasis;
using riptide::geometry::smearDirection;

constexpr std::uint64_t kProjectSeed = 20260907;
// Rounding of a handful of double products: a few ulp of 1.
constexpr double kTolerance = 1e-14;

// Axes covering both branches of the helper choice in perpendicularBasis.
const std::array<Direction, 5> kAxes{Direction{1.0, 0.0, 0.0}, Direction{0.0, 1.0, 0.0},
                                     Direction{0.0, 0.0, -1.0}, Direction{0.6, 0.0, 0.8},
                                     Direction{0.95, 0.0, std::sqrt(1.0 - (0.95 * 0.95))}};

} // namespace

TEST_CASE("perpendicularBasis gives a right-handed orthonormal frame", "[direction]") {
    for (const auto& axis : kAxes) {
        const auto [e1, e2] = perpendicularBasis(axis);
        CHECK_THAT(dot(e1, e1), WithinAbs(1.0, kTolerance));
        CHECK_THAT(dot(e2, e2), WithinAbs(1.0, kTolerance));
        CHECK_THAT(dot(e1, e2), WithinAbs(0.0, kTolerance));
        CHECK_THAT(dot(e1, axis), WithinAbs(0.0, kTolerance));
        const Direction e3 = cross(e1, e2);
        CHECK_THAT(dot(e3, axis), WithinAbs(1.0, kTolerance));
    }
}

TEST_CASE("directionAroundAxis at theta = 0 is the axis", "[direction]") {
    for (const auto& axis : kAxes) {
        const Direction d = directionAroundAxis(axis, 0.0, 1.3);
        CHECK_THAT(angleBetween(d, axis), WithinAbs(0.0, kTolerance));
    }
}

TEST_CASE("angleBetween matches known angles", "[direction]") {
    const Direction x{1.0, 0.0, 0.0};
    const Direction y{0.0, 1.0, 0.0};
    CHECK_THAT(angleBetween(x, y), WithinAbs(riptide::units::pi / 2.0, kTolerance));
    CHECK_THAT(angleBetween(x, Direction{-1.0, 0.0, 0.0}),
               WithinAbs(riptide::units::pi, kTolerance));
    for (const auto& axis : kAxes) {
        const double theta = 0.37;
        CHECK_THAT(angleBetween(axis, directionAroundAxis(axis, theta, 2.1)),
                   WithinAbs(theta, kTolerance));
    }
}

TEST_CASE("smearDirection returns unit vectors", "[direction]") {
    riptide::random::Engine engine{kProjectSeed};
    constexpr std::size_t kDraws = 1000;
    constexpr double kSigma = 0.08;
    for (std::size_t i = 0; i < kDraws; ++i) {
        const Direction d = smearDirection(kAxes[i % kAxes.size()], kSigma, engine);
        CHECK_THAT(dot(d, d), WithinAbs(1.0, kTolerance));
    }
}
