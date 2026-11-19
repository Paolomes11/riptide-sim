#pragma once

// Chord of the pencil beam through the scintillator box, and the depth of a
// vertex along the beam. The expected depth profile and the first-interaction
// fractions need the path length L that the beam travels inside the box; it
// is computed here with the generic slab method rather than assuming the beam
// runs along z, so an oblique or off-axis configuration is still checked
// correctly.

#include <optional>

#include "riptide/sim/run_config.hpp"

namespace riptide::sim::checks {

// Ray parameters (distances along the unit direction, [mm]) where the beam
// enters and leaves the box; t_exit - t_enter is the chord length L.
struct BeamChord {
    double t_enter = 0.0;
    double t_exit = 0.0;
    [[nodiscard]] double length() const {
        return t_exit - t_enter;
    }
};

// Intersection of the ray position + t * direction (t >= 0) with the
// axis-aligned box [-half_size, +half_size] centred at the origin, as the
// scintillator is built. Returns nullopt when the ray misses the box, grazes
// it with zero length, or the box lies behind the start point. Throws
// std::invalid_argument if the direction is not a unit vector (within
// kDirectionNormTolerance), a half size is not positive, or the start point
// lies strictly inside the box: the configuration forbids a gun inside the
// scintillator, and the depth would have no entry point to refer to.
[[nodiscard]] std::optional<BeamChord> rayBoxChord(const Vec3& position, const Vec3& direction,
                                                   const Vec3& half_size);

// Depth of a point along the beam, measured from the entry point:
// (point - position) . direction - t_enter [mm]. For a vertex inside the box
// it lies in [0, L].
[[nodiscard]] double depthAlongBeam(const Vec3& point, const Vec3& position, const Vec3& direction,
                                    const BeamChord& chord);

// Whether a point lies in the box enlarged by `tolerance` [mm] on every side;
// the tolerance absorbs the rounding of vertices sitting on a face.
[[nodiscard]] bool insideBox(const Vec3& point, const Vec3& half_size, double tolerance);

} // namespace riptide::sim::checks
