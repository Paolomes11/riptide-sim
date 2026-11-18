#include "riptide/sim/beam_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace riptide::sim::checks {

namespace {

// A direction component smaller than this is treated as parallel to the
// slab: dividing by it would only produce huge, rounding-dominated t values.
constexpr double kParallelTolerance = 1.0e-12;

std::array<double, 3> components(const Vec3& v) {
    return {v.x, v.y, v.z};
}

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

} // namespace

std::optional<BeamChord> rayBoxChord(const Vec3& position, const Vec3& direction,
                                     const Vec3& half_size) {
    if (std::abs(std::sqrt(dot(direction, direction)) - 1.0) > kDirectionNormTolerance) {
        throw std::invalid_argument("rayBoxChord: direction is not a unit vector");
    }
    const auto p = components(position);
    const auto d = components(direction);
    const auto h = components(half_size);
    double t_min = -std::numeric_limits<double>::infinity();
    double t_max = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < p.size(); ++i) {
        if (!(h.at(i) > 0.0)) {
            throw std::invalid_argument("rayBoxChord: half size must be positive");
        }
        if (std::abs(d.at(i)) < kParallelTolerance) {
            // Parallel to this pair of faces: either always between them or never.
            if (std::abs(p.at(i)) > h.at(i)) {
                return std::nullopt;
            }
            continue;
        }
        double t1 = (-h.at(i) - p.at(i)) / d.at(i);
        double t2 = (h.at(i) - p.at(i)) / d.at(i);
        if (t1 > t2) {
            std::swap(t1, t2);
        }
        t_min = std::max(t_min, t1);
        t_max = std::min(t_max, t2);
    }
    if (!(t_max > t_min) || t_max <= 0.0) {
        return std::nullopt;
    }
    if (t_min < 0.0) {
        throw std::invalid_argument("rayBoxChord: the beam starts inside the box");
    }
    return BeamChord{t_min, t_max};
}

double depthAlongBeam(const Vec3& point, const Vec3& position, const Vec3& direction,
                      const BeamChord& chord) {
    const Vec3 offset{point.x - position.x, point.y - position.y, point.z - position.z};
    return dot(offset, direction) - chord.t_enter;
}

bool insideBox(const Vec3& point, const Vec3& half_size, double tolerance) {
    return std::abs(point.x) <= half_size.x + tolerance &&
           std::abs(point.y) <= half_size.y + tolerance &&
           std::abs(point.z) <= half_size.z + tolerance;
}

} // namespace riptide::sim::checks
