#pragma once

// Pseudo-experiment calibration of Case C (D47): the per-experiment results
// of the two-stage analysis (D46), their ensemble numbers for mu_E, log
// sigma_E and omega_n, and the shift of the hyperparameters under the
// alternative prior uniform in sigma_E. The definitions are those of the toy
// checklist: mu_E is estimated by its posterior mean, log sigma_E by its
// posterior median, omega_n by the stage-1' MAP; the pulls use the posterior
// means. Experiments draw from their own engines and are written by index,
// so the results do not depend on the number of threads. Pure computation.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/stats/calibration_b.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_b_posterior.hpp"
#include "riptide/stats/case_c.hpp"
#include "riptide/stats/case_c_posterior.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_c.hpp"

namespace riptide::stats {

struct CaseCExperimentResult {
    CaseCTruth truth;
    CaseCHyperSummary hyper;                    // baseline prior, uniform in log sigma_E
    std::optional<CaseCHyperSummary> alternate; // prior uniform in sigma_E, if run
    PredictiveEnergy predictive;
    geometry::Direction omega_0;
    geometry::Direction omega_hat;
    DirectionSummary direction; // stage 1' against the truth
};

// analyseCaseCExperiment on the reconstructed sample, with the uniform
// direction prior and the baseline hyperparameter prior; with_alternate adds
// the summaries under the prior uniform in sigma_E on the same base.
[[nodiscard]] CaseCExperimentResult
analyseCaseCPseudoExperiment(const CaseCModel& model, const CaseCPseudoExperiment& experiment,
                             const std::vector<geometry::Direction>& coarse_grid,
                             bool with_alternate);

// One row of the checklist: `experiments` experiments of `events` events,
// experiment k on the ids of (section, k), truth drawn from `truth` and
// E_n following `spectrum`.
struct CaseCRowSpec {
    std::size_t events = 0;
    std::size_t experiments = 0;
    std::uint64_t section = 0;
    CaseCTruthRange truth = kCaseCChecklistTruth;
    CaseCSpectrum spectrum = CaseCSpectrum::truncated;
    bool with_alternate = false;
};

[[nodiscard]] std::vector<CaseCExperimentResult>
runCaseCRow(const CaseCModel& model, const EnergyDomain& domain, const CaseCRowSpec& row,
            const CaseBResolution& resolution, std::uint64_t seed,
            const std::vector<geometry::Direction>& coarse_grid, unsigned threads);

// The frequentist behaviour of one hyperparameter over a set of experiments.
// residual = estimate - truth; pull = (posterior mean - truth) / posterior
// std; pull_width is the population std of the pulls (ddof = 0), as in the
// toy.
struct HyperEnsemble {
    std::size_t entries = 0;
    double bias = 0.0;       // mean residual
    double bias_error = 0.0; // sample std of the residuals / sqrt(entries)
    double rms = 0.0;        // of the residuals
    double std_mean = 0.0;   // mean posterior std
    double pull_mean = 0.0;
    double pull_mean_error = 0.0; // sample std of the pulls / sqrt(entries)
    double pull_width = 0.0;
    std::array<double, kHyperCredibleLevels.size()> coverage{}; // equal-tailed intervals
};

struct CaseCEnsemble {
    HyperEnsemble mu;
    HyperEnsemble log_sigma;
    CaseBEnsemble direction;
    DirectionBias direction_bias;
};

// Under the baseline prior or, with `alternate`, under the alternative one
// (std::invalid_argument if an experiment lacks it; the direction is the
// same in both). Needs at least two experiments; else std::invalid_argument.
[[nodiscard]] CaseCEnsemble summariseCaseCEnsemble(const std::vector<CaseCExperimentResult>& r,
                                                   bool alternate = false);

// The direction part as Case B results, for the Case B summaries.
[[nodiscard]] std::vector<CaseBExperimentResult>
caseCDirectionResults(const std::vector<CaseCExperimentResult>& r);

// Shift (alternate - baseline) / baseline std of the estimates of mu_E (mean)
// and log sigma_E (median), as the toy robustness test: mean, its error
// (sample std / sqrt(entries)) and the largest absolute value.
struct HyperShift {
    double mean = 0.0;
    double mean_error = 0.0;
    double max_abs = 0.0;
};

struct CaseCPriorShift {
    HyperShift mu;
    HyperShift log_sigma;
};

// Needs at least two experiments, all with the alternate summaries; else
// std::invalid_argument.
[[nodiscard]] CaseCPriorShift caseCPriorShift(const std::vector<CaseCExperimentResult>& r);

} // namespace riptide::stats
