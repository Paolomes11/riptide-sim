#include "riptide/stats/validate.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>

namespace riptide::stats {

namespace {

constexpr double kHalf = 0.5;

// exp(log_post - max) normalised to 1.
std::vector<double> normalisedWeights(std::size_t n, const double* log_post) {
    const double max_value = *std::max_element(log_post, log_post + n);
    if (!std::isfinite(max_value)) {
        throw std::invalid_argument("a log-posterior needs a finite maximum");
    }
    std::vector<double> w(n);
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        w[i] = std::exp(log_post[i] - max_value);
        sum += w[i];
    }
    for (double& value : w) {
        value /= sum;
    }
    return w;
}

void requireSameSize(std::size_t a, std::size_t b) {
    if (a != b) {
        throw std::invalid_argument("validation inputs have different lengths");
    }
}

} // namespace

PosteriorMoments posteriorMeanStd(const Grid& grid, const double* log_post) {
    const std::vector<double> w = normalisedWeights(grid.size(), log_post);
    double mean = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) {
        mean += w[i] * grid[i];
    }
    double variance = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) {
        const double dx = grid[i] - mean;
        variance += w[i] * dx * dx;
    }
    return {mean, std::sqrt(variance)};
}

PosteriorMoments posteriorMeanStd(const Grid& grid, const std::vector<double>& log_post) {
    requireSameSize(grid.size(), log_post.size());
    return posteriorMeanStd(grid, log_post.data());
}

std::vector<Interval> credibleIntervals(const Grid& grid, const double* log_post,
                                        const std::vector<double>& levels) {
    const std::size_t n = grid.size();
    const std::vector<double> w = normalisedWeights(n, log_post);
    std::vector<double> cdf(n);
    double cumulative = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        cumulative += w[i];
        cdf[i] = cumulative - (kHalf * w[i]);
    }

    const auto quantile = [&](double q) {
        // First centre with cdf >= q, kept in [1, n - 1] so the bracket
        // [upper - 1, upper] always exists.
        const auto first = std::lower_bound(cdf.begin(), cdf.end(), q);
        const std::size_t upper = std::clamp(
            static_cast<std::size_t>(std::distance(cdf.begin(), first)), std::size_t{1}, n - 1);
        const double cdf_lo = cdf[upper - 1];
        const double cdf_hi = cdf[upper];
        const double span = cdf_hi - cdf_lo;
        // An empty cell (weight underflowed to 0) has no slope to invert.
        const double frac = span > 0.0 ? std::clamp((q - cdf_lo) / span, 0.0, 1.0) : 0.0;
        return grid[upper - 1] + (frac * (grid[upper] - grid[upper - 1]));
    };

    std::vector<Interval> intervals;
    intervals.reserve(levels.size());
    for (const double level : levels) {
        if (!(level > 0.0 && level < 1.0)) {
            throw std::invalid_argument("credible levels must lie in (0, 1)");
        }
        const double tail = kHalf * (1.0 - level);
        intervals.push_back({quantile(tail), quantile(1.0 - tail)});
    }
    return intervals;
}

std::vector<double> coverage(const std::vector<double>& truth, const IntervalTable& intervals) {
    requireSameSize(truth.size(), intervals.events());
    if (truth.empty()) {
        throw std::invalid_argument("coverage of an empty sample");
    }
    std::vector<double> fraction(intervals.levels.size(), 0.0);
    for (std::size_t l = 0; l < fraction.size(); ++l) {
        std::size_t covered = 0;
        for (std::size_t k = 0; k < truth.size(); ++k) {
            covered += intervals.at(k, l).contains(truth[k]) ? 1U : 0U;
        }
        fraction[l] = static_cast<double>(covered) / static_cast<double>(truth.size());
    }
    return fraction;
}

std::vector<double> pulls(const std::vector<double>& truth, const std::vector<double>& estimate,
                          const std::vector<double>& sigma_hat) {
    requireSameSize(truth.size(), estimate.size());
    requireSameSize(truth.size(), sigma_hat.size());
    std::vector<double> result(truth.size());
    for (std::size_t k = 0; k < truth.size(); ++k) {
        if (!(sigma_hat[k] > 0.0)) {
            throw std::invalid_argument("a pull needs a positive declared sigma");
        }
        result[k] = (estimate[k] - truth[k]) / sigma_hat[k];
    }
    return result;
}

BiasCurve biasCurve(const std::vector<double>& truth, const std::vector<double>& estimate,
                    std::size_t n_bins) {
    requireSameSize(truth.size(), estimate.size());
    if (truth.empty() || n_bins == 0) {
        throw std::invalid_argument("a bias curve needs events and at least one bin");
    }
    const auto [min_it, max_it] = std::minmax_element(truth.begin(), truth.end());
    const Grid edges(*min_it, *max_it, n_bins + 1);

    BiasCurve curve;
    curve.centers.resize(n_bins);
    curve.counts.assign(n_bins, 0);
    std::vector<double> sums(n_bins, 0.0);
    for (std::size_t b = 0; b < n_bins; ++b) {
        curve.centers[b] = kHalf * (edges[b] + edges[b + 1]);
    }
    for (std::size_t k = 0; k < truth.size(); ++k) {
        // numpy.digitize: the number of edges <= x, minus one, clipped.
        const auto above = std::upper_bound(edges.values().begin(), edges.values().end(), truth[k]);
        const auto index = std::distance(edges.values().begin(), above) - 1;
        const std::size_t bin =
            std::min(static_cast<std::size_t>(std::max<std::ptrdiff_t>(index, 0)), n_bins - 1);
        sums[bin] += estimate[k] - truth[k];
        ++curve.counts[bin];
    }
    curve.bias.resize(n_bins);
    for (std::size_t b = 0; b < n_bins; ++b) {
        curve.bias[b] = curve.counts[b] > 0 ? sums[b] / static_cast<double>(curve.counts[b])
                                            : std::numeric_limits<double>::quiet_NaN();
    }
    return curve;
}

Checklist runChecklist(const std::vector<double>& truth, const std::vector<double>& estimate,
                       const std::vector<double>& sigma_hat, const IntervalTable& intervals,
                       std::size_t n_bins) {
    Checklist result;
    result.bias = biasCurve(truth, estimate, n_bins);
    const std::vector<double> pull = pulls(truth, estimate, sigma_hat);
    double sum = 0.0;
    for (const double value : pull) {
        sum += value;
    }
    result.pull_mean = sum / static_cast<double>(pull.size());
    double squares = 0.0;
    for (const double value : pull) {
        const double d = value - result.pull_mean;
        squares += d * d;
    }
    result.pull_width = std::sqrt(squares / static_cast<double>(pull.size()));
    result.levels = intervals.levels;
    result.coverage = coverage(truth, intervals);
    return result;
}

} // namespace riptide::stats
