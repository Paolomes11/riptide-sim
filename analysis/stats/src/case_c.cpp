#include "riptide/stats/case_c.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "riptide/stats/log_sum_exp.hpp"
#include "riptide/stats/priors.hpp"

namespace riptide::stats {

namespace {

constexpr double kHalf = 0.5;

// Below this a dot product of the linear-space stage 2 is redone exactly.
// Each weight that underflowed to 0 was below 2^-1074 relative to its row
// maximum, and the prior weights are at most 1, so the lost part of a dot
// product is below n_En * 2^-1074; against a dot product of at least 2^-900
// that is a relative error below 1e-49, far under the rounding of the sum.
constexpr double kLinearDotFloor = 0x1p-900;

// Four independent partial sums: without fast-math the compiler keeps a
// single accumulator in order, and that chain of dependent additions would
// bound the speed; four chains run in parallel and can be vectorised. The
// order differs from a plain loop only in the last bits.
double dotProduct(const double* a, const double* b, std::size_t n) noexcept {
    double s0 = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    double s3 = 0.0;
    std::size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        s0 += a[i] * b[i];
        s1 += a[i + 1] * b[i + 1];
        s2 += a[i + 2] * b[i + 2];
        s3 += a[i + 3] * b[i + 3];
    }
    for (; i < n; ++i) {
        s0 += a[i] * b[i];
    }
    return (s0 + s1) + (s2 + s3);
}

void checkGridPoints(std::size_t n_mu, std::size_t n_sigma) {
    if (n_mu < 2 || n_sigma < 2) {
        throw std::invalid_argument("a hyperparameter grid needs at least 2 points per axis");
    }
}

} // namespace

std::vector<double> truncatedNormalEnergyLogPrior(const Grid& En_grid, double mu_E,
                                                  double sigma_E) {
    if (!std::isfinite(mu_E) || !std::isfinite(sigma_E) || !(sigma_E > 0.0)) {
        throw std::invalid_argument("the spectrum needs a finite mu_E and a positive sigma_E");
    }
    const std::size_t n = En_grid.size();
    std::vector<double> log_density(n);
    std::vector<double> weighted(n);
    const double log_step = std::log(En_grid.step());
    // Exponents relative to the node nearest mu_E,
    // (E - mu)^2 - (E_ref - mu)^2 = (E - E_ref)(E + E_ref - 2 mu): with mu_E
    // outside the domain and a small sigma_E the plain exponents reach -1e4
    // or below, and their rounding (~1e-12 absolute) would spoil the
    // normalisation; the difference form stays near 0 where the mass is.
    const double E_ref = std::clamp(mu_E, En_grid.front(), En_grid.back());
    const double inv_variance = 1.0 / (sigma_E * sigma_E);
    for (std::size_t i = 0; i < n; ++i) {
        const double E = En_grid[i];
        log_density[i] = -kHalf * (E - E_ref) * (E + E_ref - (2.0 * mu_E)) * inv_variance;
        // Trapezoid weights: the end points count half.
        const double log_weight = (i == 0 || i + 1 == n) ? log_step - std::log(2.0) : log_step;
        weighted[i] = log_density[i] + log_weight;
    }
    // The constant 1 / (sqrt(2 pi) sigma) cancels in the normalisation, so it
    // is never computed.
    const double log_norm = logSumExp(weighted);
    for (double& value : log_density) {
        value -= log_norm;
    }
    return log_density;
}

double HyperGrid::sigmaAt(std::size_t index) const noexcept {
    return std::exp(log_sigma[index % log_sigma.size()]);
}

HyperGrid globalHyperGrid(const EnergyDomain& domain, std::size_t n_mu, std::size_t n_sigma) {
    checkGridPoints(n_mu, n_sigma);
    return {Grid{domain.En_min, domain.En_max, n_mu},
            Grid{std::log(kHyperSigmaMin), std::log(kHyperSigmaMax), n_sigma}};
}

HyperWindow hyperWindow(const std::vector<double>& log_post, const HyperGrid& grid,
                        double delta_log) {
    if (log_post.size() != grid.size()) {
        throw std::invalid_argument("the log-posterior does not match the hyperparameter grid");
    }
    if (!std::isfinite(delta_log) || !(delta_log > 0.0)) {
        throw std::invalid_argument("the window threshold must be positive and finite");
    }
    const double max_value = *std::max_element(log_post.begin(), log_post.end());
    if (!std::isfinite(max_value)) {
        throw std::invalid_argument("the hyperparameter log-posterior has no finite maximum");
    }
    std::size_t i_lo = grid.mu.size();
    std::size_t i_hi = 0;
    std::size_t j_lo = grid.log_sigma.size();
    std::size_t j_hi = 0;
    for (std::size_t i = 0; i < grid.mu.size(); ++i) {
        for (std::size_t j = 0; j < grid.log_sigma.size(); ++j) {
            if (log_post[grid.index(i, j)] > max_value - delta_log) {
                i_lo = std::min(i_lo, i);
                i_hi = std::max(i_hi, i);
                j_lo = std::min(j_lo, j);
                j_hi = std::max(j_hi, j);
            }
        }
    }
    HyperWindow window;
    window.mu_lo = std::max(grid.mu[i_lo] - grid.mu.step(), grid.mu.front());
    window.mu_hi = std::min(grid.mu[i_hi] + grid.mu.step(), grid.mu.back());
    window.log_sigma_lo =
        std::max(grid.log_sigma[j_lo] - grid.log_sigma.step(), grid.log_sigma.front());
    window.log_sigma_hi =
        std::min(grid.log_sigma[j_hi] + grid.log_sigma.step(), grid.log_sigma.back());
    return window;
}

HyperGrid windowHyperGrid(const HyperWindow& window, std::size_t n_mu, std::size_t n_sigma) {
    checkGridPoints(n_mu, n_sigma);
    return {Grid{window.mu_lo, window.mu_hi, n_mu},
            Grid{window.log_sigma_lo, window.log_sigma_hi, n_sigma}};
}

std::vector<double> uniformLogSigmaHyperPrior(const HyperGrid& grid) {
    return std::vector<double>(grid.size(), -std::log(static_cast<double>(grid.size())));
}

std::vector<double> uniformSigmaHyperPrior(const HyperGrid& grid) {
    // log sigma_j - log sum_(i, j) sigma_j, summed in log space.
    const std::vector<double>& log_sigma = grid.log_sigma.values();
    const double log_norm = logSumExp(log_sigma) + std::log(static_cast<double>(grid.mu.size()));
    std::vector<double> result(grid.size());
    for (std::size_t index = 0; index < grid.size(); ++index) {
        result[index] = log_sigma[index % log_sigma.size()] - log_norm;
    }
    return result;
}

void checkHyperLogPosterior(const std::vector<double>& log_posterior) {
    bool any_finite = false;
    for (const double value : log_posterior) {
        if (std::isnan(value) || (std::isinf(value) && value > 0.0)) {
            throw std::domain_error("hyperparameter log-posterior: NaN or +inf entry");
        }
        any_finite = any_finite || std::isfinite(value);
    }
    if (!any_finite) {
        throw std::domain_error("hyperparameter log-posterior: -inf on every point, no spectrum "
                                "is compatible with all events");
    }
}

CaseCModel::CaseCModel(Grid En_grid, Grid theta_grid, Grid theta_obs_grid,
                       CaseBResolution resolution)
    : stage1_(En_grid, std::move(theta_grid), std::move(theta_obs_grid),
              flatEnergyLogPrior(En_grid), resolution),
      log_flat_prior_(flatEnergyLogPrior(En_grid).front()) {}

HierarchicalBase CaseCModel::base(const std::vector<double>& Ep_hat,
                                  const std::vector<geometry::Direction>& tracks,
                                  const geometry::Direction& omega_0) const {
    if (Ep_hat.size() != tracks.size()) {
        throw std::invalid_argument("Case C needs one track per Ep_hat");
    }
    HierarchicalBase result;
    result.events = Ep_hat.size();
    result.points = energyGrid().size();
    result.values.resize(result.events * result.points);
    for (std::size_t k = 0; k < result.events; ++k) {
        const std::vector<double> row =
            stage1_.energyLogPosterior(Ep_hat[k], geometry::angleBetween(tracks[k], omega_0));
        double* out = result.values.data() + (k * result.points);
        for (std::size_t i = 0; i < result.points; ++i) {
            out[i] = row[i] - log_flat_prior_;
        }
    }
    return result;
}

std::vector<double> CaseCModel::logLikelihood(const HierarchicalBase& base,
                                              const HyperGrid& grid) const {
    const Grid& En_grid = energyGrid();
    const std::size_t n_En = En_grid.size();
    if (base.points != n_En || base.values.size() != base.events * n_En) {
        throw std::invalid_argument("the Case C base is not on the E_n grid of the model");
    }
    // log sum_i exp(b_ki + p_i) = m_k + M + log sum_i w_ki q_i, with
    // w_ki = exp(b_ki - m_k) and q_i = exp(p_i - M) taken once per row and
    // once per point: the sum over E_n becomes a dot product instead of n_En
    // exponentials per event and point.
    std::vector<double> row_max(base.events);
    std::vector<double> weights(base.values.size());
    for (std::size_t k = 0; k < base.events; ++k) {
        const double* row = base.row(k);
        row_max[k] = *std::max_element(row, row + n_En);
        double* w = weights.data() + (k * n_En);
        for (std::size_t i = 0; i < n_En; ++i) {
            // An all -inf row keeps zero weights; it takes the exact path.
            w[i] = std::isinf(row_max[k]) ? 0.0 : std::exp(row[i] - row_max[k]);
        }
    }
    // The trapezoid weights that normalise the prior also weigh the sum, so
    // that w_i pi_i is a mass summing to 1 on every node, end nodes included.
    const double log_step = std::log(En_grid.step());
    std::vector<double> log_weight(n_En, log_step);
    log_weight.front() = log_step - std::log(2.0);
    log_weight.back() = log_step - std::log(2.0);
    std::vector<double> result(grid.size());
    std::vector<double> prior_weights(n_En);
    std::vector<double> terms(n_En);
    for (std::size_t point = 0; point < grid.size(); ++point) {
        std::vector<double> log_prior_En =
            truncatedNormalEnergyLogPrior(En_grid, grid.muAt(point), grid.sigmaAt(point));
        for (std::size_t i = 0; i < n_En; ++i) {
            log_prior_En[i] += log_weight[i];
        }
        const double prior_max = *std::max_element(log_prior_En.begin(), log_prior_En.end());
        for (std::size_t i = 0; i < n_En; ++i) {
            prior_weights[i] = std::exp(log_prior_En[i] - prior_max);
        }
        double sum = 0.0;
        for (std::size_t k = 0; k < base.events; ++k) {
            const double dot = dotProduct(weights.data() + (k * n_En), prior_weights.data(), n_En);
            if (dot >= kLinearDotFloor) {
                sum += row_max[k] + prior_max + std::log(dot);
                continue;
            }
            // The row and the prior barely overlap: the weights that
            // underflowed could carry the sum, so it is taken exactly.
            const double* row = base.row(k);
            for (std::size_t i = 0; i < n_En; ++i) {
                terms[i] = row[i] + log_prior_En[i];
            }
            // Finite or -inf by construction: the base comes from logsumexp
            // sums and the prior is finite.
            sum += logSumExpUnchecked(terms.data(), n_En);
        }
        result[point] = sum;
    }
    return result;
}

std::vector<double> CaseCModel::logPosterior(const HierarchicalBase& base, const HyperGrid& grid,
                                             const HyperLogPrior& log_prior) const {
    std::vector<double> result = logLikelihood(base, grid);
    const std::vector<double> prior = log_prior(grid);
    if (prior.size() != grid.size()) {
        throw std::invalid_argument("the hyperparameter prior does not match the grid");
    }
    for (std::size_t point = 0; point < result.size(); ++point) {
        result[point] += prior[point];
    }
    checkHyperLogPosterior(result);
    return result;
}

HyperPosterior CaseCModel::refineHyperparameters(const HierarchicalBase& base,
                                                 const HyperLogPrior& log_prior) const {
    const Grid& En_grid = energyGrid();
    const HyperGrid global =
        globalHyperGrid(EnergyDomain{En_grid.front(), En_grid.back(), En_grid.size()});
    const std::vector<double> coarse = logPosterior(base, global, log_prior);
    HyperGrid fine = windowHyperGrid(hyperWindow(coarse, global));
    std::vector<double> fine_log_post = logPosterior(base, fine, log_prior);
    return {std::move(fine), std::move(fine_log_post)};
}

} // namespace riptide::stats
