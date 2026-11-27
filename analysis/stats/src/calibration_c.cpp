#include "riptide/stats/calibration_c.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace riptide::stats {

namespace {

double mean(const std::vector<double>& v) {
    double sum = 0.0;
    for (const double x : v) {
        sum += x;
    }
    return sum / static_cast<double>(v.size());
}

// Sample std (ddof = 1) / sqrt(n), the error of the mean.
double meanError(const std::vector<double>& v, double m) {
    double sum = 0.0;
    for (const double x : v) {
        sum += (x - m) * (x - m);
    }
    const auto n = static_cast<double>(v.size());
    return std::sqrt(sum / (n - 1.0) / n);
}

// Population std (ddof = 0), the pull width of the toy checklist.
double populationStd(const std::vector<double>& v, double m) {
    double sum = 0.0;
    for (const double x : v) {
        sum += (x - m) * (x - m);
    }
    return std::sqrt(sum / static_cast<double>(v.size()));
}

double rootMeanSquare(const std::vector<double>& v) {
    double sum = 0.0;
    for (const double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum / static_cast<double>(v.size()));
}

void requireEntries(std::size_t n) {
    if (n < 2) {
        throw std::invalid_argument("Case C ensembles need at least two experiments");
    }
}

const CaseCHyperSummary& hyperOf(const CaseCExperimentResult& r, bool alternate) {
    if (!alternate) {
        return r.hyper;
    }
    if (!r.alternate) {
        throw std::invalid_argument("a Case C experiment lacks the alternate prior summaries");
    }
    return *r.alternate;
}

// estimate is the point estimate of the residual (mean for mu_E, median for
// log sigma_E); the pull always uses the posterior mean and std.
HyperEnsemble summariseHyper(const std::vector<const HyperparameterSummary*>& s,
                             const std::vector<double>& truth, bool median_estimate) {
    const std::size_t n = s.size();
    std::vector<double> residual(n);
    std::vector<double> pull(n);
    std::vector<double> std_post(n);
    std::array<std::size_t, kHyperCredibleLevels.size()> covered{};
    for (std::size_t k = 0; k < n; ++k) {
        const HyperparameterSummary& h = *s[k];
        residual[k] = (median_estimate ? h.median : h.mean) - truth[k];
        pull[k] = (h.mean - truth[k]) / h.std;
        std_post[k] = h.std;
        for (std::size_t l = 0; l < covered.size(); ++l) {
            if (h.intervals[l].contains(truth[k])) {
                ++covered[l];
            }
        }
    }
    HyperEnsemble e;
    e.entries = n;
    e.bias = mean(residual);
    e.bias_error = meanError(residual, e.bias);
    e.rms = rootMeanSquare(residual);
    e.std_mean = mean(std_post);
    e.pull_mean = mean(pull);
    e.pull_mean_error = meanError(pull, e.pull_mean);
    e.pull_width = populationStd(pull, e.pull_mean);
    for (std::size_t l = 0; l < covered.size(); ++l) {
        e.coverage[l] = static_cast<double>(covered[l]) / static_cast<double>(n);
    }
    return e;
}

HyperShift summariseShift(const std::vector<double>& shift) {
    HyperShift s;
    s.mean = mean(shift);
    s.mean_error = meanError(shift, s.mean);
    for (const double x : shift) {
        s.max_abs = std::max(s.max_abs, std::abs(x));
    }
    return s;
}

} // namespace

CaseCExperimentResult
analyseCaseCPseudoExperiment(const CaseCModel& model, const CaseCPseudoExperiment& experiment,
                             const std::vector<geometry::Direction>& coarse_grid,
                             bool with_alternate) {
    const DirectionLogPrior direction_prior = uniformDirectionLogPrior;
    const HyperLogPrior baseline = uniformLogSigmaHyperPrior;
    const HyperLogPrior alternate = uniformSigmaHyperPrior;
    const CaseBPseudoSample& sample = experiment.sample;
    CaseCAnalysis analysis =
        analyseCaseCExperiment(model, sample.reco_Ep, sample.reco_track, coarse_grid,
                               direction_prior, baseline, with_alternate ? &alternate : nullptr);
    CaseCExperimentResult r;
    r.truth = experiment.truth;
    r.hyper = analysis.hyper;
    r.alternate = analysis.alternate;
    r.predictive = analysis.predictive;
    r.omega_0 = analysis.omega_0;
    r.omega_hat = analysis.direction.omega_hat;
    r.direction = summarizeDirection(analysis.direction, experiment.truth.omega_n);
    return r;
}

std::vector<CaseCExperimentResult>
runCaseCRow(const CaseCModel& model, const EnergyDomain& domain, const CaseCRowSpec& row,
            const CaseBResolution& resolution, std::uint64_t seed,
            const std::vector<geometry::Direction>& coarse_grid, unsigned threads) {
    std::vector<CaseCExperimentResult> results(row.experiments);
    parallelForEach(row.experiments, threads, [&](std::size_t k) {
        const CaseCPseudoExperiment experiment = caseCPseudoExperiment(
            domain, row.events, resolution, row.truth, row.spectrum, seed, row.section, k);
        results[k] =
            analyseCaseCPseudoExperiment(model, experiment, coarse_grid, row.with_alternate);
    });
    return results;
}

std::vector<CaseBExperimentResult>
caseCDirectionResults(const std::vector<CaseCExperimentResult>& r) {
    std::vector<CaseBExperimentResult> out;
    out.reserve(r.size());
    for (const CaseCExperimentResult& x : r) {
        out.push_back({x.truth.omega_n, x.omega_hat, x.direction});
    }
    return out;
}

CaseCEnsemble summariseCaseCEnsemble(const std::vector<CaseCExperimentResult>& r, bool alternate) {
    requireEntries(r.size());
    const std::size_t n = r.size();
    std::vector<const HyperparameterSummary*> mu(n);
    std::vector<const HyperparameterSummary*> log_sigma(n);
    std::vector<double> mu_true(n);
    std::vector<double> log_sigma_true(n);
    for (std::size_t k = 0; k < n; ++k) {
        const CaseCHyperSummary& h = hyperOf(r[k], alternate);
        mu[k] = &h.mu;
        log_sigma[k] = &h.log_sigma;
        mu_true[k] = r[k].truth.mu;
        log_sigma_true[k] = std::log(r[k].truth.sigma);
    }
    CaseCEnsemble e;
    e.mu = summariseHyper(mu, mu_true, false);
    e.log_sigma = summariseHyper(log_sigma, log_sigma_true, true);
    const std::vector<CaseBExperimentResult> direction = caseCDirectionResults(r);
    e.direction = summariseCaseBEnsemble(direction);
    e.direction_bias = directionBias(direction);
    return e;
}

CaseCPriorShift caseCPriorShift(const std::vector<CaseCExperimentResult>& r) {
    requireEntries(r.size());
    std::vector<double> mu(r.size());
    std::vector<double> log_sigma(r.size());
    for (std::size_t k = 0; k < r.size(); ++k) {
        const CaseCHyperSummary& base = r[k].hyper;
        const CaseCHyperSummary& alt = hyperOf(r[k], true);
        mu[k] = (alt.mu.mean - base.mu.mean) / base.mu.std;
        log_sigma[k] = (alt.log_sigma.median - base.log_sigma.median) / base.log_sigma.std;
    }
    return {summariseShift(mu), summariseShift(log_sigma)};
}

} // namespace riptide::stats
