#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <stdexcept>

#include "riptide/sim/beam_geometry.hpp"

namespace checks = riptide::sim::checks;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::sim::Vec3;

namespace {

constexpr double kTolerance = 1e-12;
const Vec3 kHalfCube{50.0, 50.0, 50.0}; // the 100 mm cube of the configs

} // namespace

TEST_CASE("beam chord at normal incidence", "[beam_geometry]") {
    const Vec3 position{0.0, 0.0, -100.0};
    const Vec3 direction{0.0, 0.0, 1.0};
    const auto chord = checks::rayBoxChord(position, direction, kHalfCube);
    REQUIRE(chord.has_value());
    CHECK_THAT(chord->t_enter, WithinAbs(50.0, kTolerance));
    CHECK_THAT(chord->t_exit, WithinAbs(150.0, kTolerance));
    CHECK_THAT(chord->length(), WithinAbs(100.0, kTolerance));

    CHECK_THAT(checks::depthAlongBeam({0.0, 0.0, -50.0}, position, direction, *chord),
               WithinAbs(0.0, kTolerance));
    CHECK_THAT(checks::depthAlongBeam({0.0, 0.0, 25.0}, position, direction, *chord),
               WithinAbs(75.0, kTolerance));
}

TEST_CASE("beam chord at oblique incidence", "[beam_geometry]") {
    // Along the body diagonal from outside a corner: the chord is the full
    // diagonal of the cube, 100 sqrt(3) mm, entering after 100 sqrt(3)/2 mm.
    const double s = 1.0 / std::sqrt(3.0);
    const Vec3 direction{s, s, s};
    const Vec3 position{-100.0, -100.0, -100.0};
    const auto chord = checks::rayBoxChord(position, direction, kHalfCube);
    REQUIRE(chord.has_value());
    CHECK_THAT(chord->length(), WithinRel(100.0 * std::sqrt(3.0), kTolerance));
    CHECK_THAT(chord->t_enter, WithinRel(50.0 * std::sqrt(3.0), kTolerance));
    CHECK_THAT(checks::depthAlongBeam({0.0, 0.0, 0.0}, position, direction, *chord),
               WithinRel(50.0 * std::sqrt(3.0), kTolerance));
}

TEST_CASE("beam chord through a non-cubic box, off axis", "[beam_geometry]") {
    // A beam along x at y = 10 crosses the full x extent of the box.
    const auto chord = checks::rayBoxChord({-200.0, 10.0, 0.0}, {1.0, 0.0, 0.0}, {30.0, 20.0, 5.0});
    REQUIRE(chord.has_value());
    CHECK_THAT(chord->length(), WithinAbs(60.0, kTolerance));
    CHECK_THAT(chord->t_enter, WithinAbs(170.0, kTolerance));
}

TEST_CASE("beam that misses the box", "[beam_geometry]") {
    // Parallel to z but displaced beyond the side face.
    CHECK_FALSE(checks::rayBoxChord({60.0, 0.0, -100.0}, {0.0, 0.0, 1.0}, kHalfCube));
    // Pointing away from the box.
    CHECK_FALSE(checks::rayBoxChord({0.0, 0.0, -100.0}, {0.0, 0.0, -1.0}, kHalfCube));
    // Oblique, passing beside a corner.
    const double s = 1.0 / std::sqrt(2.0);
    CHECK_FALSE(checks::rayBoxChord({-200.0, 0.0, 0.0}, {s, s, 0.0}, kHalfCube));
}

TEST_CASE("beam chord rejects invalid input", "[beam_geometry]") {
    CHECK_THROWS_AS(checks::rayBoxChord({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, kHalfCube),
                    std::invalid_argument);
    CHECK_THROWS_AS(checks::rayBoxChord({0.0, 0.0, -100.0}, {0.0, 0.0, 2.0}, kHalfCube),
                    std::invalid_argument);
    CHECK_THROWS_AS(checks::rayBoxChord({0.0, 0.0, -100.0}, {0.0, 0.0, 1.0}, {50.0, 0.0, 50.0}),
                    std::invalid_argument);
}

TEST_CASE("point inside the box with tolerance", "[beam_geometry]") {
    CHECK(checks::insideBox({0.0, 0.0, 0.0}, kHalfCube, 0.0));
    CHECK(checks::insideBox({50.0, -50.0, 50.0}, kHalfCube, 0.0));
    CHECK_FALSE(checks::insideBox({50.001, 0.0, 0.0}, kHalfCube, 0.0));
    CHECK(checks::insideBox({50.001, 0.0, 0.0}, kHalfCube, 0.01));
}
