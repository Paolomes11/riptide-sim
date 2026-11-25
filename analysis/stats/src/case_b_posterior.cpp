#include "riptide/stats/case_b_posterior.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>

#include "riptide/common/units.hpp"
#include "riptide/stats/log_sum_exp.hpp"

namespace riptide::stats {

namespace {

// log_prior + loglik on one grid, with the prior's size checked.
std::vector<double> gridLogPosterior(const CaseBModel& model, const TrackEnergyTable& table,
                                     const std::vector<geometry::Direction>& tracks,
                                     const std::vector<geometry::Direction>& grid,
                                     const DirectionLogPrior& log_prior) {
    std::vector<double> log_post = model.summedLogLikelihood(table, tracks, grid);
    const std::vector<double> prior = log_prior(grid);
    if (prior.size() != grid.size()) {
        throw std::invalid_argument("refineSharedDirection: prior size differs from the grid");
    }
    for (std::size_t i = 0; i < log_post.size(); ++i) {
        log_post[i] += prior[i];
    }
    return log_post;
}

std::size_t argMax(const std::vector<double>& values) {
    // max_element returns the first maximum: ties resolve to the lowest index.
    return static_cast<std::size_t>(
        std::distance(values.begin(), std::max_element(values.begin(), values.end())));
}

// exp(log_post - max), normalised to unit sum.
std::vector<double> normalisedWeights(const std::vector<double>& log_post) {
    const double max_value = *std::max_element(log_post.begin(), log_post.end());
    std::vector<double> weights(log_post.size());
    double sum = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) {
        weights[i] = std::exp(log_post[i] - max_value);
        sum += weights[i];
    }
    for (double& w : weights) {
        w /= sum;
    }
    return weights;
}

} // namespace

std::vector<double> uniformDirectionLogPrior(const std::vector<geometry::Direction>& grid) {
    if (grid.empty()) {
        throw std::invalid_argument("uniformDirectionLogPrior: empty grid");
    }
    return std::vector<double>(grid.size(), -std::log(static_cast<double>(grid.size())));
}

std::vector<double> vmfDirectionLogPrior(const std::vector<geometry::Direction>& grid,
                                         const geometry::Direction& axis, double kappa) {
    if (grid.empty()) {
        throw std::invalid_argument("vmfDirectionLogPrior: empty grid");
    }
    if (!std::isfinite(kappa) || kappa < 0.0) {
        throw std::invalid_argument("vmfDirectionLogPrior: kappa must be finite and >= 0");
    }
    std::vector<double> log_w(grid.size());
    for (std::size_t i = 0; i < grid.size(); ++i) {
        log_w[i] = kappa * (geometry::dot(grid[i], axis) - 1.0);
    }
    const double log_norm = logSumExp(log_w);
    for (double& value : log_w) {
        value -= log_norm;
    }
    return log_w;
}

void checkDirectionLogPosterior(const std::vector<double>& log_posterior) {
    bool any_finite = false;
    for (const double value : log_posterior) {
        if (std::isnan(value) || (std::isinf(value) && value > 0.0)) {
            throw std::domain_error("direction log-posterior: NaN or +inf entry");
        }
        any_finite = any_finite || std::isfinite(value);
    }
    if (!any_finite) {
        throw std::domain_error("direction log-posterior: -inf on every candidate, no direction "
                                "is kinematically compatible with all tracks");
    }
}

DirectionPosterior refineSharedDirection(const CaseBModel& model, const TrackEnergyTable& table,
                                         const std::vector<geometry::Direction>& tracks,
                                         const std::vector<geometry::Direction>& coarse_grid,
                                         const DirectionLogPrior& log_prior, double delta_log,
                                         std::size_t n_cap) {
    const std::vector<double> coarse =
        gridLogPosterior(model, table, tracks, coarse_grid, log_prior);
    checkDirectionLogPosterior(coarse);
    const geometry::Direction coarse_map = coarse_grid[argMax(coarse)];

    DirectionPosterior result;
    result.cap_radius = directionCapRadius(coarse, coarse_grid, delta_log);
    result.cap = directionCapGrid(coarse_map, result.cap_radius, n_cap);
    result.cap_log_posterior = gridLogPosterior(model, table, tracks, result.cap, log_prior);
    checkDirectionLogPosterior(result.cap_log_posterior);
    result.map_index = argMax(result.cap_log_posterior);
    result.omega_hat = result.cap[result.map_index];
    return result;
}

DirectionSummary summarizeDirection(const DirectionPosterior& posterior,
                                    const geometry::Direction& omega_true) {
    const std::vector<geometry::Direction>& cap = posterior.cap;
    if (cap.empty() || cap.size() != posterior.cap_log_posterior.size()) {
        throw std::invalid_argument("summarizeDirection: empty cap or sizes that differ");
    }
    const std::vector<double> weights = normalisedWeights(posterior.cap_log_posterior);

    DirectionSummary summary;
    summary.error = geometry::angleBetween(omega_true, posterior.omega_hat);

    double second_moment = 0.0;
    double min_cos = 1.0;
    std::size_t nearest = 0;
    double nearest_cos = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < cap.size(); ++i) {
        const double angle = geometry::angleBetween(posterior.omega_hat, cap[i]);
        second_moment += weights[i] * angle * angle;
        min_cos = std::min(min_cos, geometry::dot(cap[i], posterior.omega_hat));
        const double cos_truth = geometry::dot(cap[i], omega_true);
        if (cos_truth > nearest_cos) {
            nearest_cos = cos_truth;
            nearest = i;
        }
    }
    summary.sigma = std::sqrt(second_moment);
    summary.pull = summary.error / summary.sigma;

    const double cap_pixel =
        std::sqrt(2.0 * units::pi * (1.0 - min_cos) / static_cast<double>(cap.size()));
    summary.truth_in_cap =
        geometry::angleBetween(cap[nearest], omega_true) < kTruthInCapPixels * cap_pixel;

    const double p_nearest = weights[nearest];
    double mass_above = 0.0;
    for (const double w : weights) {
        if (w > p_nearest) {
            mass_above += w;
        }
    }
    for (std::size_t l = 0; l < kDirectionCredibleLevels.size(); ++l) {
        summary.covered[l] = summary.truth_in_cap && mass_above < kDirectionCredibleLevels[l];
    }
    return summary;
}

} // namespace riptide::stats
