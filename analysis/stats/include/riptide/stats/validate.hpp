#pragma once

// Validation on simulated data with known truth, in the order of the
// checklist: bias -> resolution -> pull -> coverage. pull = (estimate - truth)
// / declared sigma is ~N(0, 1) for a calibrated reconstruction; coverage is
// the fraction of level-L credible intervals that contain the truth and must
// be ~L. Contraction with N and prior robustness are run by the calibration
// app on top of these functions.

#include <cstddef>
#include <vector>

#include "riptide/stats/grid.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

struct PosteriorMoments {
    double mean = 0.0;
    double std = 0.0;
};

// Mean and standard deviation of an unnormalised log-posterior on a grid,
// with weights exp(log_post - max) normalised to 1. The log-posterior must
// have one finite maximum, else std::invalid_argument.
[[nodiscard]] PosteriorMoments posteriorMeanStd(const Grid& grid, const double* log_post);
[[nodiscard]] PosteriorMoments posteriorMeanStd(const Grid& grid,
                                                const std::vector<double>& log_post);

struct Interval {
    double lo = 0.0;
    double hi = 0.0;

    // End points included, as in the toy coverage.
    [[nodiscard]] bool contains(double value) const noexcept {
        return lo <= value && value <= hi;
    }
};

// Equal-tailed credible intervals, one per level in (0, 1). The CDF is taken
// at the cell centres (cumulative mass minus half the cell) and inverted by
// linear interpolation between grid points; quantiles beyond the first or
// last centre stop at the grid edge.
[[nodiscard]] std::vector<Interval> credibleIntervals(const Grid& grid, const double* log_post,
                                                      const std::vector<double>& levels);

// The 0.5 quantile, with the same CDF as credibleIntervals: the limit of its
// interval as the level goes to 0.
[[nodiscard]] double posteriorMedian(const Grid& grid, const double* log_post);

// Credible intervals of many events: entry [k * levels.size() + l] is event
// k at level l.
struct IntervalTable {
    std::vector<double> levels;
    std::vector<Interval> values;

    [[nodiscard]] std::size_t events() const noexcept {
        return levels.empty() ? 0 : values.size() / levels.size();
    }
    [[nodiscard]] const Interval& at(std::size_t k, std::size_t l) const noexcept {
        return values[(k * levels.size()) + l];
    }
};

// Fraction of the intervals at each level that contain the truth.
[[nodiscard]] std::vector<double> coverage(const std::vector<double>& truth,
                                           const IntervalTable& intervals);

// (estimate - truth) / sigma_hat, event by event.
[[nodiscard]] std::vector<double> pulls(const std::vector<double>& truth,
                                        const std::vector<double>& estimate,
                                        const std::vector<double>& sigma_hat);

struct BiasCurve {
    std::vector<double> centers;
    std::vector<double> bias; // NaN in an empty bin
    std::vector<std::size_t> counts;
};

// Mean residual estimate - truth in n_bins equal bins spanning the range of
// the truth. Bin edges follow numpy.linspace; an event on an inner edge goes
// to the bin above it and the maximum goes to the last bin.
[[nodiscard]] BiasCurve biasCurve(const std::vector<double>& truth,
                                  const std::vector<double>& estimate, std::size_t n_bins);

inline constexpr std::size_t kDefaultBiasBins = 20;

struct Checklist {
    BiasCurve bias;
    double pull_mean = 0.0;
    double pull_width = 0.0; // population standard deviation of the pulls
    std::vector<double> levels;
    std::vector<double> coverage;
};

// The first four checklist items on one set of reconstructed events.
[[nodiscard]] Checklist runChecklist(const std::vector<double>& truth,
                                     const std::vector<double>& estimate,
                                     const std::vector<double>& sigma_hat,
                                     const IntervalTable& intervals,
                                     std::size_t n_bins = kDefaultBiasBins);

} // namespace riptide::stats
