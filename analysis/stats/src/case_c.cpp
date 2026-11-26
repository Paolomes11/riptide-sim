#include "riptide/stats/case_c.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "riptide/stats/log_sum_exp.hpp"

namespace riptide::stats {

namespace {

constexpr double kHalf = 0.5;

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

} // namespace riptide::stats
