#pragma once

// Posterior of the shared neutron direction of Case B and its per-experiment
// summaries (D42). The posterior is log prior + sum_k loglik_k on the
// equal-area coarse sphere, then again on a cap of the same lattice around
// the coarse MAP; the estimate is the MAP on the cap. The prior is counted
// once per grid, never once per event.

#include <array>
#include <cstddef>
#include <functional>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// Credible levels of the HPD coverage, as in the toy model.
inline constexpr std::array<double, 3> kDirectionCredibleLevels{0.68, 0.90, 0.95};
// The truth counts as inside the cap only within this many cap pixels of its
// nearest pixel: farther away it lies outside the region where the posterior
// was computed and cannot be covered.
inline constexpr double kTruthInCapPixels = 3.0;

// Log prior weights of a set of equal-area pixels; the weights must sum to one
// on that set. Called once on the coarse grid and once on the cap.
using DirectionLogPrior =
    std::function<std::vector<double>(const std::vector<geometry::Direction>&)>;

// -log n on every pixel: the uniform prior in solid angle, since all pixels
// have the same area. Throws std::invalid_argument on an empty grid.
[[nodiscard]] std::vector<double>
uniformDirectionLogPrior(const std::vector<geometry::Direction>& grid);

// von Mises-Fisher prior of unit axis `axis` and concentration kappa >= 0,
// normalised over the pixels: log w_i = kappa (cos alpha_i - 1) - log sum.
// The shift by -1 keeps exp() from overflowing at large kappa; kappa = 0 is
// the uniform prior. Only for the robustness test of the prior. Throws
// std::invalid_argument on an empty grid or a negative or non-finite kappa.
[[nodiscard]] std::vector<double> vmfDirectionLogPrior(const std::vector<geometry::Direction>& grid,
                                                       const geometry::Direction& axis,
                                                       double kappa);

struct DirectionPosterior {
    geometry::Direction omega_hat;         // MAP on the cap
    std::size_t map_index = 0;             // index of omega_hat in cap
    double cap_radius = 0.0;               // rad
    std::vector<geometry::Direction> cap;  // refinement pixels
    std::vector<double> cap_log_posterior; // unnormalised, one per pixel
};

// Throws std::domain_error if a log-posterior has a NaN or +inf entry or no
// finite entry: an all -inf posterior would otherwise give pixel 0 as MAP
// without warning (no direction is kinematically compatible with every
// track).
void checkDirectionLogPosterior(const std::vector<double>& log_posterior);

// Coarse posterior on `coarse_grid`, cap of radius directionCapRadius about
// its MAP with n_cap pixels, posterior on the cap with the same prior
// function, MAP on the cap. Ties resolve to the lowest index. Throws as
// checkDirectionLogPosterior on either grid, and std::invalid_argument for
// mismatched table and tracks (summedLogLikelihood) or a prior of the wrong
// size.
[[nodiscard]] DirectionPosterior
refineSharedDirection(const CaseBModel& model, const TrackEnergyTable& table,
                      const std::vector<geometry::Direction>& tracks,
                      const std::vector<geometry::Direction>& coarse_grid,
                      const DirectionLogPrior& log_prior, double delta_log = kCapDeltaLog,
                      std::size_t n_cap = kCapGridPixels);

struct DirectionSummary {
    double error = 0.0; // angle(omega_true, omega_hat), rad
    double sigma = 0.0; // sqrt(E_post[angle(omega_hat, Omega)^2]), rad
    double pull = 0.0;  // error / sigma
    bool truth_in_cap = false;
    std::array<bool, kDirectionCredibleLevels.size()> covered{};
};

// Summaries of one experiment, as the toy checklist computes them. The
// posterior weights are exp(log_post - max) normalised on the cap. The pixel
// size is sqrt(2 pi (1 - min_i cos(cap_i, omega_hat)) / n_cap). Coverage at
// level L is HPD on the equal-area cap: the mass of the pixels strictly more
// probable than the one nearest the truth is below L, and the truth is in the
// cap.
[[nodiscard]] DirectionSummary summarizeDirection(const DirectionPosterior& posterior,
                                                  const geometry::Direction& omega_true);

} // namespace riptide::stats
