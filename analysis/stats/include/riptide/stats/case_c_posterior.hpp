#pragma once

// Summaries of the Case C posterior and the full two-stage analysis of one
// experiment (D46): stage 1 with the flat E_n prior gives omega_0, stage 2
// gives the spectrum hyperparameters on their fine grid, and stage 1' fits the
// direction again with the posterior predictive spectrum as E_n prior. The
// hyperparameters and the direction are never put on one joint grid.

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/stats/case_b_posterior.hpp"
#include "riptide/stats/case_c.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/validate.hpp"

namespace riptide::stats {

// Credible levels of the equal-tailed hyperparameter intervals, the same as
// those of the direction so that one table holds all coverages.
inline constexpr std::array<double, 3> kHyperCredibleLevels = kDirectionCredibleLevels;

// Normalised log marginals on the fine grid: of mu_E (over its axis, summing
// out log sigma_E) and of log sigma_E (summing out mu_E). Throws
// std::invalid_argument if the log-posterior does not match the grid.
[[nodiscard]] std::vector<double> muLogMarginal(const HyperPosterior& posterior);
[[nodiscard]] std::vector<double> logSigmaLogMarginal(const HyperPosterior& posterior);

// One hyperparameter: posterior mean, standard deviation and median, and the
// equal-tailed intervals at kHyperCredibleLevels.
struct HyperparameterSummary {
    double mean = 0.0;
    double std = 0.0;
    double median = 0.0;
    std::array<Interval, kHyperCredibleLevels.size()> intervals{};
};

// mu_E in MeV, log sigma_E in log(MeV).
struct CaseCHyperSummary {
    HyperparameterSummary mu;
    HyperparameterSummary log_sigma;
};

[[nodiscard]] HyperparameterSummary
summarizeHyperparameter(const Grid& axis, const std::vector<double>& log_marginal);
[[nodiscard]] CaseCHyperSummary summarizeHyperparameters(const HyperPosterior& posterior);

// Mean and standard deviation of E_n under the posterior predictive
// E_n ~ N(mu_E, sigma_E) averaged over p(mu_E, sigma_E | D):
//   E[E_n] = E[mu_E],  Var[E_n] = E[sigma_E^2] + Var[mu_E].
// A plug-in of the point estimates would ignore the uncertainty on the
// spectrum and make the stage-1' direction over-confident at small N.
struct PredictiveEnergy {
    double mean = 0.0; // MeV
    double std = 0.0;  // MeV
};

[[nodiscard]] PredictiveEnergy predictiveEnergyMoments(const HyperPosterior& posterior);

// Stage 1': a Case B model on the grids and with the acceptance of stage 1,
// with the truncated N(mean, std) prior on E_n, then refineSharedDirection. Throws as
// truncatedNormalEnergyLogPrior and refineSharedDirection.
[[nodiscard]] DirectionPosterior
refineHierarchicalDirection(const CaseCModel& model, const std::vector<double>& Ep_hat,
                            const std::vector<geometry::Direction>& tracks,
                            const std::vector<geometry::Direction>& coarse_grid,
                            const DirectionLogPrior& log_prior, const PredictiveEnergy& spectrum);

struct CaseCAnalysis {
    geometry::Direction omega_0;                // stage-1 MAP, flat E_n prior
    CaseCHyperSummary hyper;                    // baseline hyperparameter prior
    std::optional<CaseCHyperSummary> alternate; // robustness prior, if given
    PredictiveEnergy predictive;                // from the baseline posterior
    DirectionPosterior direction;               // stage 1'
};

// The chain of D46 on one experiment. The predictive spectrum of stage 1'
// always comes from the baseline prior; the alternate prior only adds its
// hyperparameter summaries, so the two priors are compared on the same base.
[[nodiscard]] CaseCAnalysis
analyseCaseCExperiment(const CaseCModel& model, const std::vector<double>& Ep_hat,
                       const std::vector<geometry::Direction>& tracks,
                       const std::vector<geometry::Direction>& coarse_grid,
                       const DirectionLogPrior& direction_prior, const HyperLogPrior& hyper_prior,
                       const HyperLogPrior* alternate_prior = nullptr);

} // namespace riptide::stats
