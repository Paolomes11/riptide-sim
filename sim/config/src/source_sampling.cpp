#include "riptide/sim/source_sampling.hpp"

#include <algorithm>
#include <cmath>

namespace riptide::sim {

namespace {

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Faces of the box: each pair has area (2 h_a)(2 h_b); summed with the
// weight |cos| of its normal, one face of each pair faces the beam.
constexpr double kFaceAreaFactor = 4.0;
// Surface of a box: 2 (ab + ac + bc) with full edges, 8 (h_a h_b + ...) with
// half edges; the mean projected area of a convex body is a quarter of it.
constexpr double kSurfaceFactor = 8.0;
constexpr double kCauchyMeanShadow = 0.25;

} // namespace

OrthonormalPair orthonormalPair(const Vec3& n) {
    const double sign = std::copysign(1.0, n.z);
    const double a = -1.0 / (sign + n.z);
    const double b = n.x * n.y * a;
    return {{1.0 + sign * n.x * n.x * a, sign * b, -sign * n.x}, {b, sign + n.y * n.y * a, -n.y}};
}

Vec3 isotropicDirection(double u1, double u2) {
    const double cos_theta = 1.0 - 2.0 * u1;
    // max() guards the rounding of 1 - cos^2 below zero at the poles.
    const double sin_theta = std::sqrt(std::max(0.0, (1.0 - cos_theta) * (1.0 + cos_theta)));
    const double phi = kTwoPi * u2;
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
}

double halfDiagonal(const Vec3& half_size) {
    return std::sqrt(dot(half_size, half_size));
}

double projectedBoxArea(const Vec3& direction, const Vec3& half_size) {
    const auto& h = half_size;
    return kFaceAreaFactor *
           (std::abs(direction.x) * h.y * h.z + std::abs(direction.y) * h.x * h.z +
            std::abs(direction.z) * h.x * h.y);
}

double shadowRadius(const Vec3& direction, const Vec3& half_size) {
    // The distance of a corner c from the line is sqrt(|c|^2 - (c . d)^2);
    // |c| is the same for all corners, so the farthest corner is the one with
    // the smallest |c . d|, which the signs of the corner choose independently.
    const auto& h = half_size;
    double min_projection = halfDiagonal(h);
    for (const double sy : {-1.0, 1.0}) {
        for (const double sz : {-1.0, 1.0}) {
            const Vec3 corner{h.x, sy * h.y, sz * h.z};
            min_projection = std::min(min_projection, std::abs(dot(corner, direction)));
        }
    }
    const double diagonal = halfDiagonal(h);
    return std::sqrt(std::max(0.0, (diagonal - min_projection) * (diagonal + min_projection)));
}

FarFieldDisk farFieldDisk(const SourceConfig& source, const Vec3& half_size) {
    if (source.beam != BeamMode::FarField) {
        throw ConfigError("farFieldDisk: the source is not far-field");
    }
    const double diagonal = halfDiagonal(half_size);
    const double radius = source.isotropic ? diagonal : shadowRadius(source.direction, half_size);
    return {radius, diagonal + source.margin};
}

Vec3 farFieldStart(const Vec3& direction, const FarFieldDisk& disk, double u1, double u2) {
    const auto basis = orthonormalPair(direction);
    const double rho = disk.radius * std::sqrt(u1);
    const double phi = kTwoPi * u2;
    const double a = rho * std::cos(phi);
    const double b = rho * std::sin(phi);
    return {-disk.distance * direction.x + a * basis.e1.x + b * basis.e2.x,
            -disk.distance * direction.y + a * basis.e1.y + b * basis.e2.y,
            -disk.distance * direction.z + a * basis.e1.z + b * basis.e2.z};
}

double farFieldAcceptance(const SourceConfig& source, const Vec3& half_size) {
    const auto disk = farFieldDisk(source, half_size);
    const double disk_area = units::pi * disk.radius * disk.radius;
    const auto& h = half_size;
    const double shadow =
        source.isotropic ? kCauchyMeanShadow * kSurfaceFactor * (h.x * h.y + h.x * h.z + h.y * h.z)
                         : projectedBoxArea(source.direction, half_size);
    return shadow / disk_area;
}

double sourceExtent(const SourceConfig& source, const Vec3& half_size) {
    if (source.beam == BeamMode::Pencil) {
        return std::max({std::abs(source.position.x), std::abs(source.position.y),
                         std::abs(source.position.z)});
    }
    const auto disk = farFieldDisk(source, half_size);
    return std::hypot(disk.distance, disk.radius);
}

} // namespace riptide::sim
