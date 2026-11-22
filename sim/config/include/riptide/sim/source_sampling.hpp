#pragma once

// Geometry and sampling of the neutron source (D27). Pure functions with no
// Geant4 dependency: the generator calls them with Geant4's engine, the tests
// with a seeded standard engine.
//
// The far-field source is a parallel beam from a distant source: for each
// event a direction Omega (fixed, or uniform on the sphere) and a start point
// uniform on a disk perpendicular to Omega, centred on the line through the
// box centre and placed upstream at the half-diagonal plus a margin, so no
// start point lies inside the box. Every box point then sees the same
// direction and the same fluence, which is the far-field assumption of the
// analysis.

#include <cmath>

#include "riptide/common/units.hpp"
#include "riptide/sim/run_config.hpp"

namespace riptide::sim {

inline constexpr double kTwoPi = 2.0 * units::pi;

// Two unit vectors that complete the unit vector n to a right-handed
// orthonormal basis (e1, e2, n). Branch-free construction of Duff et al.,
// "Building an Orthonormal Basis, Revisited", JCGT 6(1), 2017, which stays
// accurate for every n, including n close to -z.
struct OrthonormalPair {
    Vec3 e1;
    Vec3 e2;
};
[[nodiscard]] OrthonormalPair orthonormalPair(const Vec3& n);

// Unit vector from two uniforms in [0, 1): cos(theta) = 1 - 2 u1 is uniform
// in (-1, 1] and phi = 2 pi u2, which is uniform on the sphere.
[[nodiscard]] Vec3 isotropicDirection(double u1, double u2);

// Half-diagonal of the box: the radius of its circumscribed sphere [mm].
[[nodiscard]] double halfDiagonal(const Vec3& half_size);

// Area of the shadow of the box on a plane perpendicular to the unit vector
// direction [mm^2]: each pair of opposite faces contributes the area of one
// face times |cos| of its normal with the direction.
[[nodiscard]] double projectedBoxArea(const Vec3& direction, const Vec3& half_size);

// Largest distance of a box corner from the line through the centre along
// the unit vector direction: the smallest centred disk that holds the shadow
// of the box [mm].
[[nodiscard]] double shadowRadius(const Vec3& direction, const Vec3& half_size);

// Start disk of the far-field source: radius and distance of its centre from
// the box centre, upstream along -direction [mm].
struct FarFieldDisk {
    double radius = 0.0;
    double distance = 0.0;
};

// For a fixed direction the disk is the smallest one that holds the shadow;
// for an isotropic source it is the half-diagonal for every direction, so the
// directions of the neutrons that hit the box keep the weight of their
// projected area, as for a true isotropic flux (D27). Throws ConfigError if
// the source is not far-field.
[[nodiscard]] FarFieldDisk farFieldDisk(const SourceConfig& source, const Vec3& half_size);

// Start point on the disk from two uniforms in [0, 1): radius r sqrt(u1) and
// angle 2 pi u2, uniform in area.
[[nodiscard]] Vec3 farFieldStart(const Vec3& direction, const FarFieldDisk& disk, double u1,
                                 double u2);

// Expected fraction of far-field neutrons whose line crosses the box: the
// projected area over the disk area, averaged over the directions for an
// isotropic source (mean projected area of a convex body = surface / 4,
// Cauchy).
[[nodiscard]] double farFieldAcceptance(const SourceConfig& source, const Vec3& half_size);

// Upper bound of |start point| over all events of the source [mm].
[[nodiscard]] double sourceExtent(const SourceConfig& source, const Vec3& half_size);

// Direction and start point of one far-field event. `uniform` returns doubles
// in [0, 1); it is called twice for the direction when the source is
// isotropic, then twice for the start point, in this order.
struct SourceRay {
    Vec3 position;
    Vec3 direction;
};
template <class Uniform>
[[nodiscard]] SourceRay sampleFarField(const SourceConfig& source, const FarFieldDisk& disk,
                                       Uniform&& uniform) {
    SourceRay ray;
    if (source.isotropic) {
        const double u1 = uniform();
        const double u2 = uniform();
        ray.direction = isotropicDirection(u1, u2);
    } else {
        ray.direction = source.direction;
    }
    const double u3 = uniform();
    const double u4 = uniform();
    ray.position = farFieldStart(ray.direction, disk, u3, u4);
    return ray;
}

} // namespace riptide::sim
