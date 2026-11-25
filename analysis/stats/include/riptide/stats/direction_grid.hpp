#pragma once

// Grids on the unit sphere for the shared neutron direction of Case B (D42).
// The coarse grid covers the whole sphere with a Fibonacci lattice; the MAP is
// then refined on a cap of the same lattice around the coarse MAP. Both
// lattices give every pixel the same solid angle, so a uniform prior per pixel
// is uniform in solid angle.

#include <cstddef>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// Defaults of the toy model, so that its numbers compare at the same
// discretisation.
inline constexpr std::size_t kSphereGridPixels = 3000;
inline constexpr std::size_t kCapGridPixels = 8000;
// Coarse pixels whose log-posterior is within this of the maximum are kept
// inside the cap: exp(-10) ~ 4.5e-5 of the peak density.
inline constexpr double kCapDeltaLog = 10.0;

// n pixels of equal area: cos theta_i = 1 - 2 (i + 1/2) / n, phi_i = i times
// the golden angle pi (3 - sqrt 5). Throws std::invalid_argument for n == 0.
[[nodiscard]] std::vector<geometry::Direction> sphereGrid(std::size_t n);

// n pixels of equal area on the cap of half-opening `radius` about the unit
// vector `centre`: cos alpha_i = 1 - (1 - cos radius) (i + 1/2) / n, azimuth
// i times the golden angle. Throws std::invalid_argument for n == 0 or a
// radius outside (0, pi].
[[nodiscard]] std::vector<geometry::Direction> directionCapGrid(const geometry::Direction& centre,
                                                                double radius, std::size_t n);

// Radius of the refinement cap: the largest angle from the coarse MAP of a
// pixel with log-posterior above max - delta_log, plus two coarse steps (the
// angle from the MAP to its nearest pixel), at most pi. The MAP is excluded
// from the step by index, not by angle > 0: the rounded angle of a pixel to
// itself is ~1e-8 rad, not 0, and would give a degenerate cap. Throws
// std::invalid_argument for sizes that differ, fewer than two pixels, a
// non-finite maximum, a NaN entry or a non-positive delta_log.
[[nodiscard]] double directionCapRadius(const std::vector<double>& log_posterior,
                                        const std::vector<geometry::Direction>& grid,
                                        double delta_log = kCapDeltaLog);

} // namespace riptide::stats
