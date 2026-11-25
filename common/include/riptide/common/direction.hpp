#pragma once

// Unit vectors in 3D and the frame operations on them shared by the detector
// response (analysis/io) and the direction posteriors (analysis/stats). They
// live here, header-only, because analysis/stats may not depend on io.
//
// A direction is a plain (x, y, z) triple; the functions below expect unit
// vectors where stated and do not renormalise their output beyond what the
// construction guarantees.

#include <cmath>
#include <utility>

#include "riptide/common/random.hpp"

namespace riptide::geometry {

struct Direction {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

namespace detail {

// Above this |axis.x| the helper axis switches from x to y, so that the cross
// product with the helper never comes close to zero.
inline constexpr double kHelperSwitch = 0.9;

} // namespace detail

[[nodiscard]] inline Direction cross(const Direction& a, const Direction& b) noexcept {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

[[nodiscard]] inline double dot(const Direction& a, const Direction& b) noexcept {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]] inline Direction normalized(const Direction& v) noexcept {
    const double norm = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return {v.x / norm, v.y / norm, v.z / norm};
}

// Angle between two unit vectors in [0, pi]. atan2(|a x b|, a.b) keeps full
// precision near 0 and pi, where acos(a.b) loses half the digits.
[[nodiscard]] inline double angleBetween(const Direction& a, const Direction& b) noexcept {
    const Direction c = cross(a, b);
    return std::atan2(std::sqrt(dot(c, c)), dot(a, b));
}

// Two unit vectors (e1, e2) such that (e1, e2, axis) is a right-handed
// orthonormal frame; `axis` must be a unit vector.
[[nodiscard]] inline std::pair<Direction, Direction>
perpendicularBasis(const Direction& axis) noexcept {
    const Direction helper = std::abs(axis.x) < detail::kHelperSwitch ? Direction{1.0, 0.0, 0.0}
                                                                      : Direction{0.0, 1.0, 0.0};
    const Direction e1 = normalized(cross(axis, helper));
    return {e1, cross(axis, e1)};
}

// Unit vector at polar angle theta and azimuth phi about the unit vector `axis`.
[[nodiscard]] inline Direction directionAroundAxis(const Direction& axis, double theta,
                                                   double phi) noexcept {
    const auto [e1, e2] = perpendicularBasis(axis);
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);
    const double along_e1 = sin_theta * std::cos(phi);
    const double along_e2 = sin_theta * std::sin(phi);
    return {cos_theta * axis.x + along_e1 * e1.x + along_e2 * e2.x,
            cos_theta * axis.y + along_e1 * e1.y + along_e2 * e2.y,
            cos_theta * axis.z + along_e1 * e1.z + along_e2 * e2.z};
}

// The track moved by a Gaussian offset in the tangent plane, two draws of
// N(0, sigma_theta); `track` must be a unit vector, the result is one.
[[nodiscard]] inline Direction smearDirection(const Direction& track, double sigma_theta,
                                              random::Engine& engine) {
    // Two draws in a fixed order (first e1, then e2): with a function call per
    // argument the order of evaluation would be unspecified.
    const double offset_1 = sigma_theta * random::standardNormal(engine);
    const double offset_2 = sigma_theta * random::standardNormal(engine);
    return directionAroundAxis(track, std::hypot(offset_1, offset_2),
                               std::atan2(offset_2, offset_1));
}

} // namespace riptide::geometry
