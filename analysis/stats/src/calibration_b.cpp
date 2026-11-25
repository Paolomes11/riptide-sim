#include "riptide/stats/calibration_b.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>

#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/calibration.hpp"
#include "riptide/stats/direction_grid.hpp"

namespace riptide::stats {

namespace {

// The first failure of one thread's stripe.
struct Failure {
    std::size_t index = std::numeric_limits<std::size_t>::max();
    std::exception_ptr error;
};

void runStripe(std::size_t first, std::size_t count, std::size_t stride,
               const std::function<void(std::size_t)>& task, Failure& failure) {
    for (std::size_t k = first; k < count; k += stride) {
        try {
            task(k);
        } catch (...) {
            // Later indices of the stripe cannot lower the failing index.
            failure = {k, std::current_exception()};
            return;
        }
    }
}

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

double rootMeanSquare(const std::vector<double>& v) {
    double sum = 0.0;
    for (const double x : v) {
        sum += x * x;
    }
    return std::sqrt(sum / static_cast<double>(v.size()));
}

// As numpy.median: the mean of the two middle values for an even count.
double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const std::size_t half = v.size() / 2;
    if (v.size() % 2 == 1) {
        return v[half];
    }
    return 0.5 * (v[half - 1] + v[half]);
}

void requireEntries(std::size_t n) {
    if (n < 2) {
        throw std::invalid_argument("Case B ensembles need at least two experiments");
    }
}

} // namespace

void parallelForEach(std::size_t count, unsigned threads,
                     const std::function<void(std::size_t)>& task) {
    if (threads == 0) {
        throw std::invalid_argument("parallelForEach needs at least one thread");
    }
    const std::size_t workers = std::min<std::size_t>(threads, std::max<std::size_t>(count, 1));
    std::vector<Failure> failures(workers);
    if (workers == 1) {
        runStripe(0, count, 1, task, failures[0]);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(workers);
        for (std::size_t t = 0; t < workers; ++t) {
            pool.emplace_back(runStripe, t, count, workers, std::cref(task), std::ref(failures[t]));
        }
        for (std::thread& thread : pool) {
            thread.join();
        }
    }
    const auto first =
        std::min_element(failures.begin(), failures.end(),
                         [](const Failure& a, const Failure& b) { return a.index < b.index; });
    if (first->error) {
        std::rethrow_exception(first->error);
    }
}

CaseBExperimentResult analyseCaseBExperiment(const CaseBModel& model,
                                             const CaseBPseudoExperiment& experiment,
                                             const std::vector<geometry::Direction>& coarse_grid,
                                             const DirectionLogPrior& log_prior) {
    const CaseBPseudoSample& sample = experiment.sample;
    const TrackEnergyTable table = model.trackEnergyTable(sample.reco_Ep);
    const DirectionPosterior posterior =
        refineSharedDirection(model, table, sample.reco_track, coarse_grid, log_prior);
    return {experiment.truth_omega_n, posterior.omega_hat,
            summarizeDirection(posterior, experiment.truth_omega_n)};
}

std::vector<CaseBExperimentResult>
runCaseBRow(const CaseBModel& model, const EnergyDomain& truth_domain, const CaseBRowSpec& row,
            const CaseBResolution& resolution, std::uint64_t seed,
            const std::vector<geometry::Direction>& coarse_grid, unsigned threads) {
    std::vector<CaseBExperimentResult> results(row.experiments);
    const DirectionLogPrior uniform = uniformDirectionLogPrior;
    parallelForEach(row.experiments, threads, [&](std::size_t k) {
        const CaseBPseudoExperiment experiment =
            caseBPseudoExperiment(truth_domain, row.events, resolution, seed, row.section, k);
        results[k] = analyseCaseBExperiment(model, experiment, coarse_grid, uniform);
    });
    return results;
}

CaseBPriorRow runCaseBPriorRow(const CaseBModel& model, const CaseBModel& energy_model,
                               const EnergyDomain& truth_domain, const CaseBRowSpec& row,
                               const CaseBPriorSpec& prior, const CaseBResolution& resolution,
                               std::uint64_t seed,
                               const std::vector<geometry::Direction>& coarse_grid,
                               unsigned threads) {
    if (!(prior.kappa >= 0.0) || !std::isfinite(prior.kappa)) {
        throw std::invalid_argument("the vMF prior needs a finite kappa >= 0");
    }
    if (!(prior.axis_angle >= 0.0 && prior.axis_angle <= units::pi)) {
        throw std::invalid_argument("the vMF axis angle must lie in [0, pi]");
    }
    CaseBPriorRow out;
    out.vmf.resize(row.experiments);
    out.energy.resize(row.experiments);
    const DirectionLogPrior uniform = uniformDirectionLogPrior;
    // The experiments are regenerated from the ids of the checklist row, so
    // they are the very ones of its reference results.
    parallelForEach(row.experiments, threads, [&](std::size_t k) {
        const CaseBPseudoExperiment experiment =
            caseBPseudoExperiment(truth_domain, row.events, resolution, seed, row.section, k);
        auto engine = random::perEventEngine(seed, experimentFirstEventId(prior.axis_section, k));
        const double phi = 2.0 * units::pi * random::uniform53(engine);
        const geometry::Direction axis =
            geometry::directionAroundAxis(experiment.truth_omega_n, prior.axis_angle, phi);
        const DirectionLogPrior vmf = [&axis, &prior](const std::vector<geometry::Direction>& g) {
            return vmfDirectionLogPrior(g, axis, prior.kappa);
        };
        out.vmf[k] = analyseCaseBExperiment(model, experiment, coarse_grid, vmf);
        out.energy[k] = analyseCaseBExperiment(energy_model, experiment, coarse_grid, uniform);
    });
    return out;
}

CaseBEnsemble summariseCaseBEnsemble(const std::vector<CaseBExperimentResult>& r) {
    requireEntries(r.size());
    const std::size_t n = r.size();
    std::vector<double> error(n);
    std::vector<double> sigma(n);
    std::vector<double> pull(n);
    std::vector<double> scaled(n);
    CaseBEnsemble e;
    e.entries = n;
    std::size_t tail = 0;
    std::size_t in_cap = 0;
    std::array<std::size_t, kDirectionCredibleLevels.size()> covered{};
    for (std::size_t k = 0; k < n; ++k) {
        const DirectionSummary& s = r[k].summary;
        error[k] = s.error;
        sigma[k] = s.sigma;
        pull[k] = s.pull;
        scaled[k] = std::sqrt(2.0) * s.pull;
        tail += scaled[k] > kScaledPullTail ? 1U : 0U;
        in_cap += s.truth_in_cap ? 1U : 0U;
        for (std::size_t l = 0; l < covered.size(); ++l) {
            covered[l] += s.covered[l] ? 1U : 0U;
        }
    }
    const auto count = static_cast<double>(n);
    e.error_mean = mean(error);
    e.error_mean_error = meanError(error, e.error_mean);
    e.error_rms = rootMeanSquare(error);
    e.sigma_mean = mean(sigma);
    e.sigma_mean_error = meanError(sigma, e.sigma_mean);
    e.pull_rms = rootMeanSquare(pull);
    e.scaled_pull_median = median(scaled);
    e.scaled_pull_tail_fraction = static_cast<double>(tail) / count;
    e.truth_in_cap_fraction = static_cast<double>(in_cap) / count;
    for (std::size_t l = 0; l < covered.size(); ++l) {
        e.coverage[l] = static_cast<double>(covered[l]) / count;
    }
    return e;
}

double logLogSlope(const std::vector<double>& events, const std::vector<double>& rms) {
    if (events.size() != rms.size() || events.size() < 2) {
        throw std::invalid_argument("logLogSlope needs two or more matching points");
    }
    std::vector<double> x(events.size());
    std::vector<double> y(rms.size());
    for (std::size_t k = 0; k < x.size(); ++k) {
        if (!(events[k] > 0.0 && rms[k] > 0.0)) {
            throw std::invalid_argument("logLogSlope needs positive values");
        }
        x[k] = std::log(events[k]);
        y[k] = std::log(rms[k]);
    }
    const double mx = mean(x);
    const double my = mean(y);
    double sxx = 0.0;
    double sxy = 0.0;
    for (std::size_t k = 0; k < x.size(); ++k) {
        sxx += (x[k] - mx) * (x[k] - mx);
        sxy += (x[k] - mx) * (y[k] - my);
    }
    if (!(sxx > 0.0)) {
        throw std::invalid_argument("logLogSlope needs at least two distinct N");
    }
    return sxy / sxx;
}

std::array<double, 2> tangentCoordinates(const geometry::Direction& estimate,
                                         const geometry::Direction& truth) {
    const double c = geometry::dot(estimate, truth);
    const geometry::Direction perp{estimate.x - (c * truth.x), estimate.y - (c * truth.y),
                                   estimate.z - (c * truth.z)};
    const double norm = std::sqrt(geometry::dot(perp, perp));
    // At norm 0 the estimate is the truth (or its antipode, where the
    // direction of the map is undefined): the coordinates are left at 0.
    const double scale = norm > 0.0 ? geometry::angleBetween(estimate, truth) / norm : 0.0;
    const double theta = std::atan2(std::hypot(truth.x, truth.y), truth.z);
    const double phi = std::atan2(truth.y, truth.x);
    const geometry::Direction e_theta{std::cos(theta) * std::cos(phi),
                                      std::cos(theta) * std::sin(phi), -std::sin(theta)};
    const geometry::Direction e_phi{-std::sin(phi), std::cos(phi), 0.0};
    return {scale * geometry::dot(perp, e_theta), scale * geometry::dot(perp, e_phi)};
}

DirectionPriorShift directionPriorShift(const std::vector<CaseBExperimentResult>& reference,
                                        const std::vector<CaseBExperimentResult>& alternative) {
    if (reference.size() != alternative.size()) {
        throw std::invalid_argument("prior shift: the two sets differ in size");
    }
    requireEntries(reference.size());
    const std::size_t n = reference.size();
    std::vector<double> shift(n);
    std::vector<double> sigma(n);
    std::array<std::vector<double>, 2> d{std::vector<double>(n), std::vector<double>(n)};
    for (std::size_t k = 0; k < n; ++k) {
        const CaseBExperimentResult& ref = reference[k];
        const CaseBExperimentResult& alt = alternative[k];
        if (!(geometry::angleBetween(ref.omega_true, alt.omega_true) <= kSameTruthAngle)) {
            throw std::invalid_argument("prior shift: the two sets have different truths");
        }
        shift[k] = geometry::angleBetween(ref.omega_hat, alt.omega_hat) / ref.summary.sigma;
        sigma[k] = ref.summary.sigma;
        const auto t_ref = tangentCoordinates(ref.omega_hat, ref.omega_true);
        const auto t_alt = tangentCoordinates(alt.omega_hat, ref.omega_true);
        d[0][k] = t_alt[0] - t_ref[0];
        d[1][k] = t_alt[1] - t_ref[1];
    }
    DirectionPriorShift out;
    out.shift_mean = mean(shift);
    out.shift_error = meanError(shift, out.shift_mean);
    out.shift_max = *std::max_element(shift.begin(), shift.end());
    const double scale = mean(sigma);
    for (std::size_t a = 0; a < 2; ++a) {
        const double m = mean(d[a]);
        out.vector_mean[a] = m / scale;
        out.vector_error[a] = meanError(d[a], m) / scale;
    }
    return out;
}

} // namespace riptide::stats
