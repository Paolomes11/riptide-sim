#pragma once

// Angular resolution of the proton track in Case B (D41). The observed track
// is a von Mises-Fisher draw about the true one, with concentration
// kappa = 1 / sigma_theta^2: to O(sigma_theta^4) the same as the Gaussian
// offset in the tangent plane of geometry::smearDirection. The likelihood
// depends on the true track only through its angle theta to the neutron
// direction, so the kernel below is the vMF density already integrated over
// the azimuth of the true track about that direction.

#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// Above this argument log I0e uses the asymptotic series, below it the
// standard library Bessel function. Both are within ~2e-15 of the reference
// values at the switch, and the series is the more accurate of the two above
// it (the library function overflows past x ~ 700).
inline constexpr double kBesselAsymptoticFrom = 30.0;

// log(I0(x) exp(-x)) for x >= 0, +inf included (-> -inf). For x above
// kBesselAsymptoticFrom: I0(x) ~ e^x / sqrt(2 pi x) sum_k c_k / (8x)^k with
// c_k = prod_{j<=k} (2j - 1)^2 / k! (Abramowitz and Stegun 9.7.1 at nu = 0).
// Throws std::invalid_argument for a negative or NaN x.
[[nodiscard]] double logI0e(double x);

// kappa = 1 / sigma_theta^2. Throws std::invalid_argument unless sigma_theta
// is finite and positive.
[[nodiscard]] double vmfConcentration(double sigma_theta);

// log K(theta_obs | theta), the density per steradian of the observed track
// at angle theta_obs from the neutron direction, per unit of the true angle
// theta, integrated over the azimuth of the true track:
//
//   K = kappa / (1 - e^{-2 kappa}) sin(theta) exp(-2 kappa sin^2((theta_obs - theta) / 2))
//       I0e(kappa sin(theta_obs) sin(theta)).
//
// For any theta, int_0^pi sin(theta_obs) K dtheta_obs = sin(theta). Tracks
// past the pole are inside I0, so no mirrored term is needed. -inf at
// theta = 0, where the azimuthal circle has no length. Angles in [0, pi];
// no checks, as it fills tables in loops (kappa from vmfConcentration).
[[nodiscard]] double logMeridianVmfKernel(double theta_obs, double theta, double kappa) noexcept;

// log(cos(theta) / pi): the non-relativistic density per steradian of the
// recoil track about the neutron direction for isotropic CM scattering. The
// E_n-dependent relativistic factor is kept apart, in the energy term (D41).
// -inf for theta > pi/2; at the rounded pi/2 of a grid, cos is ~6e-17 and the
// value is finite but negligible, as in the toy model.
[[nodiscard]] double logTrackDensity(double theta) noexcept;

} // namespace riptide::stats
