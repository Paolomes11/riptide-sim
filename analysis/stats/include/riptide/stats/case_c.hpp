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

#include "riptide/common/direction.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/case_b.hpp"
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

// Throws std::domain_error if a hyperparameter log-posterior has a NaN or
// +inf entry or no finite entry: an all -inf posterior would otherwise give
// point 0 as MAP without warning.
void checkHyperLogPosterior(const std::vector<double>& log_posterior);

// Per-event log-likelihood on the E_n grid with omega_n fixed, theta_p
// marginalised and no E_n prior, row-major: row k is event k. It does not
// depend on the hyperparameters, so it is computed once per experiment and
// reused on the global and the fine grid.
struct HierarchicalBase {
    std::size_t events = 0;
    std::size_t points = 0;
    std::vector<double> values;

    [[nodiscard]] const double* row(std::size_t k) const noexcept {
        return values.data() + (k * points);
    }
};

// The refined stage 2: the fine grid on the window and the unnormalised
// log-posterior on its points.
struct HyperPosterior {
    HyperGrid grid;
    std::vector<double> log_posterior;
};

// Stages 1 and 2 of D46. Stage 1 is a Case B model with the flat E_n prior;
// its MAP direction omega_0 fixes the angles theta_obs of stage 2, where
//
//   log L(mu_E, sigma_E) = sum_k log sum_i w_i exp[base_k(E_n,i) + log pi(E_n,i | mu_E, sigma_E)]
//
// with base_k the Case B energy posterior of event k at theta_obs minus the
// flat prior and w_i the trapezoid weights, the same rule that normalises pi.
// Without w_i the end nodes would count twice their mass: a sigma_E below the
// E_n step puts all of pi on one node, and on an end node each event would
// gain log 2, which at N = 300 pulls the posterior of a spectrum near the
// domain edge onto (mu_E at the edge, sigma_E at its minimum).
class CaseCModel {
public:
    // The grids and the resolution of the stage-1 Case B model; throws as
    // the CaseBModel constructor. log_acceptance is empty for the ideal
    // detector, or the caseBAcceptanceLogWeights table of D40 on the same
    // grids. That table is q(theta_p | E_n, sel), normalised in theta_p at
    // each E_n, so it leaves the E_n marginal alone: the spectrum of stage 2
    // stays the detected one, the density the hyperparameters describe.
    CaseCModel(Grid En_grid, Grid theta_grid, Grid theta_obs_grid, CaseBResolution resolution,
               const std::vector<double>& log_acceptance = {});

    [[nodiscard]] const CaseBModel& stage1Model() const noexcept {
        return stage1_;
    }
    [[nodiscard]] const Grid& energyGrid() const noexcept {
        return stage1_.energyGrid();
    }

    // base_k on the E_n grid for every event, with theta_obs the angle between
    // track k and omega_0. Throws std::invalid_argument for mismatched Ep_hat
    // and tracks or a non-finite Ep_hat. Unit vectors are expected.
    [[nodiscard]] HierarchicalBase base(const std::vector<double>& Ep_hat,
                                        const std::vector<geometry::Direction>& tracks,
                                        const geometry::Direction& omega_0) const;

    // log L at every point of `grid`, prior not included. Throws
    // std::invalid_argument if the base is not on the E_n grid.
    [[nodiscard]] std::vector<double> logLikelihood(const HierarchicalBase& base,
                                                    const HyperGrid& grid) const;

    // log prior + log L, the prior counted once. Throws std::invalid_argument
    // for a prior of the wrong size and as checkHyperLogPosterior.
    [[nodiscard]] std::vector<double> logPosterior(const HierarchicalBase& base,
                                                   const HyperGrid& grid,
                                                   const HyperLogPrior& log_prior) const;

    // logPosterior on the global grid, then on the fine grid of its window;
    // the prior function is evaluated on each grid.
    [[nodiscard]] HyperPosterior refineHyperparameters(const HierarchicalBase& base,
                                                       const HyperLogPrior& log_prior) const;

private:
    CaseBModel stage1_;
    // The flat log prior inside the stage-1 energy posterior, removed from
    // the base.
    double log_flat_prior_;
};

} // namespace riptide::stats
