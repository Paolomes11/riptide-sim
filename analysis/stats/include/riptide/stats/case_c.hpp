#pragma once

// Case C (D46): a shared direction omega_n and per-event energies drawn from
// one Gaussian spectrum N(mu_E, sigma_E) truncated to the E_n domain. The
// hyperparameters are inferred on a 2D grid with omega_n fixed by a Case B
// fit first: a joint 4D grid (2 direction + 2 spectrum parameters) is never
// built, because the direction and the spectrum are nearly decoupled and a 4D
// grid at the needed resolution does not fit in memory or time.

#include <cstddef>
#include <functional>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// Global hyperparameter grid: mu_E linear on the E_n domain, sigma_E
// geometric because it is a scale parameter. 0.01 MeV is far below
// sigma_Ep = 0.1 MeV, so the grid reaches the limit where the spread between
// events cannot be told from the measurement noise (Case A); 50 MeV is five
// times the domain width, where the truncated Gaussian is flat on the domain
// (Case B).
inline constexpr std::size_t kHyperMuPoints = 60;
inline constexpr std::size_t kHyperSigmaPoints = 60;
inline constexpr double kHyperSigmaMin = 0.01 * units::MeV;
inline constexpr double kHyperSigmaMax = 50.0 * units::MeV;

// The fine grid covers the window of points within kHyperWindowDeltaLog of
// the global maximum: exp(-10) ~ 5e-5 of the peak density, so the mass left
// out is negligible for moments and 95% intervals.
inline constexpr std::size_t kHyperFinePoints = 40;
inline constexpr double kHyperWindowDeltaLog = 10.0;

// log pi(E_n | mu_E, sigma_E) on the E_n grid: the Gaussian normalised by its
// trapezoidal integral over the grid, i.e. truncated to the domain, as a
// density in 1/MeV. Computed in log space, so a small sigma_E far from mu_E
// gives large negative values instead of an underflow to log(0); mu_E may lie
// outside the domain. Throws std::invalid_argument unless mu_E is finite and
// sigma_E is positive and finite.
[[nodiscard]] std::vector<double> truncatedNormalEnergyLogPrior(const Grid& En_grid, double mu_E,
                                                                double sigma_E);

// Product grid on (mu_E, log sigma_E), flattened with mu_E outer: point
// index(i, j) is (mu[i], exp(log_sigma[j])).
struct HyperGrid {
    Grid mu;
    Grid log_sigma;

    [[nodiscard]] std::size_t size() const noexcept {
        return mu.size() * log_sigma.size();
    }
    [[nodiscard]] std::size_t index(std::size_t i_mu, std::size_t j_sigma) const noexcept {
        return (i_mu * log_sigma.size()) + j_sigma;
    }
    [[nodiscard]] double muAt(std::size_t index) const noexcept {
        return mu[index / log_sigma.size()];
    }
    [[nodiscard]] double sigmaAt(std::size_t index) const noexcept;
};

// The global grid: mu_E over the E_n domain, sigma_E over
// [kHyperSigmaMin, kHyperSigmaMax].
[[nodiscard]] HyperGrid globalHyperGrid(const EnergyDomain& domain,
                                        std::size_t n_mu = kHyperMuPoints,
                                        std::size_t n_sigma = kHyperSigmaPoints);

// Bounds of a hyperparameter window, mu_E in MeV and sigma_E in log(MeV).
struct HyperWindow {
    double mu_lo = 0.0;
    double mu_hi = 0.0;
    double log_sigma_lo = 0.0;
    double log_sigma_hi = 0.0;
};

// The box holding every point of `grid` with log_post > max - delta_log,
// widened by one grid step on each side (a log step for sigma_E) so that the
// fine grid reaches past the last kept point, then clipped to the box of
// `grid`. Throws std::invalid_argument unless log_post matches the grid and
// has a finite maximum.
[[nodiscard]] HyperWindow hyperWindow(const std::vector<double>& log_post, const HyperGrid& grid,
                                      double delta_log = kHyperWindowDeltaLog);

// n_mu x n_sigma points spanning the window, built like the global grid.
[[nodiscard]] HyperGrid windowHyperGrid(const HyperWindow& window,
                                        std::size_t n_mu = kHyperFinePoints,
                                        std::size_t n_sigma = kHyperFinePoints);

// Proper prior on the points of a hyperparameter grid, as log masses summing
// to 1. Passed to the model as an argument: the prior is data, never a
// branch inside the model.
using HyperLogPrior = std::function<std::vector<double>(const HyperGrid&)>;

// The same mass at every point. On the geometric sigma_E axis this is the
// prior uniform in (mu_E, log sigma_E), the standard choice for a location
// and a scale parameter (the baseline of D46).
[[nodiscard]] std::vector<double> uniformLogSigmaHyperPrior(const HyperGrid& grid);

// Mass proportional to sigma_E: on the geometric axis a cell is ~sigma_E wide,
// so this is the prior uniform in (mu_E, sigma_E). Only for the prior
// robustness check.
[[nodiscard]] std::vector<double> uniformSigmaHyperPrior(const HyperGrid& grid);

} // namespace riptide::stats
