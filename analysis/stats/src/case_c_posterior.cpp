#include "riptide/stats/case_c_posterior.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "riptide/stats/case_b.hpp"
#include "riptide/stats/log_sum_exp.hpp"

namespace riptide::stats {

namespace {

void requireMatchingPosterior(const HyperPosterior& posterior) {
    if (posterior.log_posterior.size() != posterior.grid.size()) {
        throw std::invalid_argument("the hyperparameter log-posterior does not match its grid");
    }
}

// log sum over the other axis, then normalised: a log marginal whose
// exponentials sum to 1. `along_mu` picks the axis that is kept.
std::vector<double> logMarginal(const HyperPosterior& posterior, bool along_mu) {
    requireMatchingPosterior(posterior);
    const HyperGrid& grid = posterior.grid;
    const std::size_t n_keep = along_mu ? grid.mu.size() : grid.log_sigma.size();
    const std::size_t n_sum = along_mu ? grid.log_sigma.size() : grid.mu.size();
    std::vector<double> result(n_keep);
    std::vector<double> terms(n_sum);
    for (std::size_t a = 0; a < n_keep; ++a) {
        for (std::size_t b = 0; b < n_sum; ++b) {
            terms[b] = posterior.log_posterior[along_mu ? grid.index(a, b) : grid.index(b, a)];
        }
        result[a] = logSumExpUnchecked(terms.data(), n_sum);
    }
    const double log_norm = logSumExp(result);
    for (double& value : result) {
        value -= log_norm;
    }
    return result;
}

} // namespace

std::vector<double> muLogMarginal(const HyperPosterior& posterior) {
    return logMarginal(posterior, true);
}

std::vector<double> logSigmaLogMarginal(const HyperPosterior& posterior) {
    return logMarginal(posterior, false);
}

HyperparameterSummary summarizeHyperparameter(const Grid& axis,
                                              const std::vector<double>& log_marginal) {
    const PosteriorMoments moments = posteriorMeanStd(axis, log_marginal);
    const std::vector<Interval> intervals = credibleIntervals(
        axis, log_marginal.data(),
        std::vector<double>(kHyperCredibleLevels.begin(), kHyperCredibleLevels.end()));
    HyperparameterSummary summary;
    summary.mean = moments.mean;
    summary.std = moments.std;
    summary.median = posteriorMedian(axis, log_marginal.data());
    std::copy(intervals.begin(), intervals.end(), summary.intervals.begin());
    return summary;
}

CaseCHyperSummary summarizeHyperparameters(const HyperPosterior& posterior) {
    return {summarizeHyperparameter(posterior.grid.mu, muLogMarginal(posterior)),
            summarizeHyperparameter(posterior.grid.log_sigma, logSigmaLogMarginal(posterior))};
}

PredictiveEnergy predictiveEnergyMoments(const HyperPosterior& posterior) {
    requireMatchingPosterior(posterior);
    const std::vector<double>& log_post = posterior.log_posterior;
    const double max_value = *std::max_element(log_post.begin(), log_post.end());
    if (!std::isfinite(max_value)) {
        throw std::invalid_argument("the hyperparameter log-posterior has no finite maximum");
    }
    double total = 0.0;
    double sum_mu = 0.0;
    double sum_second = 0.0;
    for (std::size_t point = 0; point < log_post.size(); ++point) {
        const double w = std::exp(log_post[point] - max_value);
        const double mu = posterior.grid.muAt(point);
        const double sigma = posterior.grid.sigmaAt(point);
        total += w;
        sum_mu += w * mu;
        sum_second += w * ((sigma * sigma) + (mu * mu));
    }
    const double mean = sum_mu / total;
    // E[sigma^2 + mu^2] - E[mu]^2 >= E[sigma^2] > 0 up to rounding, which
    // cannot reach sigma^2 >= kHyperSigmaMin^2 = 1e-4 MeV^2.
    return {mean, std::sqrt((sum_second / total) - (mean * mean))};
}

DirectionPosterior refineHierarchicalDirection(const CaseCModel& model,
                                               const std::vector<double>& Ep_hat,
                                               const std::vector<geometry::Direction>& tracks,
                                               const std::vector<geometry::Direction>& coarse_grid,
                                               const DirectionLogPrior& log_prior,
                                               const PredictiveEnergy& spectrum) {
    const CaseBModel& stage1 = model.stage1Model();
    const CaseBModel stage1_prime(
        stage1.energyGrid(), stage1.recoilAngleGrid(), stage1.observedAngleGrid(),
        truncatedNormalEnergyLogPrior(stage1.energyGrid(), spectrum.mean, spectrum.std),
        stage1.resolution(), stage1.logAcceptance());
    const TrackEnergyTable table = stage1_prime.trackEnergyTable(Ep_hat);
    return refineSharedDirection(stage1_prime, table, tracks, coarse_grid, log_prior);
}

CaseCAnalysis analyseCaseCExperiment(const CaseCModel& model, const std::vector<double>& Ep_hat,
                                     const std::vector<geometry::Direction>& tracks,
                                     const std::vector<geometry::Direction>& coarse_grid,
                                     const DirectionLogPrior& direction_prior,
                                     const HyperLogPrior& hyper_prior,
                                     const HyperLogPrior* alternate_prior) {
    const CaseBModel& stage1 = model.stage1Model();
    const TrackEnergyTable table = stage1.trackEnergyTable(Ep_hat);
    const DirectionPosterior first =
        refineSharedDirection(stage1, table, tracks, coarse_grid, direction_prior);

    const HierarchicalBase base = model.base(Ep_hat, tracks, first.omega_hat);
    const HyperPosterior hyper = model.refineHyperparameters(base, hyper_prior);

    CaseCAnalysis result;
    result.omega_0 = first.omega_hat;
    result.hyper = summarizeHyperparameters(hyper);
    if (alternate_prior != nullptr) {
        result.alternate =
            summarizeHyperparameters(model.refineHyperparameters(base, *alternate_prior));
    }
    result.predictive = predictiveEnergyMoments(hyper);
    result.direction = refineHierarchicalDirection(model, Ep_hat, tracks, coarse_grid,
                                                   direction_prior, result.predictive);
    return result;
}

} // namespace riptide::stats
