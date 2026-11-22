#include "riptide/stats/calibration.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "riptide/common/random.hpp"
#include "riptide/stats/combine.hpp"

namespace riptide::stats {

namespace {

// Events per likelihood block in appendSingleEventSummaries: the matrix of a
// block is 1024 x n_En doubles (4 MB on the default grid), whatever the
// sample size.
constexpr std::size_t kBlockEvents = 1024;

void requireCompatible(const CaseAModel& model, const std::vector<double>& log_prior_En,
                       const PosteriorSummaries& out) {
    if (log_prior_En.size() != model.energyGrid().size()) {
        throw std::invalid_argument("the E_n prior does not match the E_n grid");
    }
    if (out.intervals.levels.empty()) {
        throw std::invalid_argument("summaries need their credible levels (makeSummaries)");
    }
}

void append(const Grid& grid, const double* log_post, double truth, PosteriorSummaries& out) {
    const PosteriorMoments moments = posteriorMeanStd(grid, log_post);
    out.truth.push_back(truth);
    out.mean.push_back(moments.mean);
    out.std.push_back(moments.std);
    for (const Interval& interval : credibleIntervals(grid, log_post, out.intervals.levels)) {
        out.intervals.values.push_back(interval);
    }
}

double mean(const std::vector<double>& values) {
    double sum = 0.0;
    for (const double v : values) {
        sum += v;
    }
    return sum / static_cast<double>(values.size());
}

// Sum of squared deviations from `centre`.
double squares(const std::vector<double>& values, double centre) {
    double sum = 0.0;
    for (const double v : values) {
        const double d = v - centre;
        sum += d * d;
    }
    return sum;
}

} // namespace

PosteriorSummaries makeSummaries(const std::vector<double>& levels) {
    if (levels.empty()) {
        throw std::invalid_argument("summaries need at least one credible level");
    }
    for (const double level : levels) {
        if (!(level > 0.0 && level < 1.0)) {
            throw std::invalid_argument("credible levels must lie in (0, 1)");
        }
    }
    PosteriorSummaries summaries;
    summaries.intervals.levels = levels;
    return summaries;
}

void appendSingleEventSummaries(const CaseAModel& model, const std::vector<double>& log_prior_En,
                                const PseudoSample& sample, PosteriorSummaries& out) {
    requireCompatible(model, log_prior_En, out);
    const Grid& grid = model.energyGrid();
    std::vector<double> log_post(grid.size());
    for (std::size_t first = 0; first < sample.size(); first += kBlockEvents) {
        const std::size_t last = std::min(first + kBlockEvents, sample.size());
        const auto lo = static_cast<std::ptrdiff_t>(first);
        const auto hi = static_cast<std::ptrdiff_t>(last);
        const std::vector<double> Ep(sample.reco_Ep.begin() + lo, sample.reco_Ep.begin() + hi);
        const std::vector<double> theta(sample.reco_theta_p.begin() + lo,
                                        sample.reco_theta_p.begin() + hi);
        const LogLikelihoodMatrix loglik = model.logLikelihoods(Ep, theta);
        for (std::size_t k = 0; k < loglik.events; ++k) {
            const double* row = loglik.row(k);
            for (std::size_t i = 0; i < grid.size(); ++i) {
                log_post[i] = row[i] + log_prior_En[i];
            }
            append(grid, log_post.data(), sample.truth_En[first + k], out);
        }
    }
}

void appendCombinedSummary(const CaseAModel& model, const std::vector<double>& log_prior_En,
                           const PseudoSample& sample, double truth_En, PosteriorSummaries& out) {
    requireCompatible(model, log_prior_En, out);
    if (sample.size() == 0) {
        throw std::invalid_argument("a combined posterior needs at least one event");
    }
    const LogLikelihoodMatrix loglik = model.logLikelihoods(sample.reco_Ep, sample.reco_theta_p);
    const std::vector<double> log_post = combineLogLikelihoods(loglik, log_prior_En);
    append(model.energyGrid(), log_post.data(), truth_En, out);
}

EnsembleSummary summariseEnsemble(const PosteriorSummaries& summaries) {
    const std::size_t n = summaries.size();
    if (n < 2) {
        throw std::invalid_argument("an ensemble summary needs at least two entries");
    }
    std::vector<double> error(n);
    for (std::size_t k = 0; k < n; ++k) {
        error[k] = summaries.mean[k] - summaries.truth[k];
    }
    const std::vector<double> pull = pulls(summaries.truth, summaries.mean, summaries.std);

    EnsembleSummary result;
    result.entries = n;
    const auto count = static_cast<double>(n);
    result.bias_mean = mean(error);
    result.bias_error = std::sqrt(squares(error, result.bias_mean) / (count - 1.0) / count);
    result.rms_error = std::sqrt(squares(error, 0.0) / count);
    result.mean_std = mean(summaries.std);
    result.std_error = std::sqrt(squares(summaries.std, result.mean_std) / (count - 1.0) / count);
    result.pull_mean = mean(pull);
    result.pull_width = std::sqrt(squares(pull, result.pull_mean) / count);
    result.levels = summaries.intervals.levels;
    result.coverage = coverage(summaries.truth, summaries.intervals);
    return result;
}

PriorShift priorShift(const PosteriorSummaries& reference, const PosteriorSummaries& alternative) {
    const std::size_t n = reference.size();
    if (n == 0 || alternative.size() != n) {
        throw std::invalid_argument("a prior shift needs two non-empty summaries of equal size");
    }
    PriorShift result;
    double shift_sum = 0.0;
    double ratio_sum = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        // Comparing floats exactly is deliberate: the two summaries must come
        // from the very same truths, not from close ones.
        if (std::isless(reference.truth[k], alternative.truth[k]) ||
            std::isgreater(reference.truth[k], alternative.truth[k])) {
            throw std::invalid_argument("a prior shift needs the same truths in both summaries");
        }
        if (!(reference.std[k] > 0.0)) {
            throw std::invalid_argument("a prior shift needs a positive reference std");
        }
        const double shift = (alternative.mean[k] - reference.mean[k]) / reference.std[k];
        shift_sum += shift;
        result.shift_max_abs = std::max(result.shift_max_abs, std::abs(shift));
        ratio_sum += alternative.std[k] / reference.std[k];
    }
    result.shift_mean = shift_sum / static_cast<double>(n);
    result.std_ratio_mean = ratio_sum / static_cast<double>(n);
    return result;
}

std::uint64_t experimentFirstEventId(std::uint64_t section, std::uint64_t index) {
    constexpr unsigned kSectionShift = kExperimentIdBits + kExperimentIndexBits;
    constexpr unsigned kIdBits = 64;
    if (index >> kExperimentIndexBits != 0 || section >> (kIdBits - kSectionShift) != 0) {
        throw std::invalid_argument("pseudo-experiment section or index out of range");
    }
    return (section << kSectionShift) | (index << kExperimentIdBits);
}

PseudoExperiment pseudoExperimentFromDomain(const EnergyDomain& truth_domain, std::size_t n,
                                            const CaseAResolution& resolution, std::uint64_t seed,
                                            std::uint64_t section, std::uint64_t index) {
    if (!(truth_domain.En_min > 0.0 && truth_domain.En_max > truth_domain.En_min)) {
        throw std::invalid_argument("pseudo-experiments need 0 < En_min < En_max");
    }
    if (n > kMaxExperimentEvents) {
        throw std::invalid_argument("too many events for one pseudo-experiment block");
    }
    const std::uint64_t first = experimentFirstEventId(section, index);
    auto engine = random::perEventEngine(seed, first + kMaxExperimentEvents);
    const double En = truth_domain.En_min +
                      ((truth_domain.En_max - truth_domain.En_min) * random::uniform53(engine));
    return {En, pseudoSampleAtEnergy(En, n, resolution, seed, first)};
}

PseudoExperiment pseudoExperimentAtEnergy(double En, std::size_t n,
                                          const CaseAResolution& resolution, std::uint64_t seed,
                                          std::uint64_t section, std::uint64_t index) {
    if (n > kMaxExperimentEvents) {
        throw std::invalid_argument("too many events for one pseudo-experiment block");
    }
    return {En,
            pseudoSampleAtEnergy(En, n, resolution, seed, experimentFirstEventId(section, index))};
}

} // namespace riptide::stats
