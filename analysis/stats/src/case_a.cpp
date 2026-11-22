#include "riptide/stats/case_a.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/log_sum_exp.hpp"

namespace riptide::stats {

namespace {

constexpr double kHalf = 0.5;

// Cut C of the theta_p window (D35 addendum). Every skipped term is at least C
// below the row maximum, so the skipped mass is at most n_theta * exp(-C)
// relative to the kept one: 500 * exp(-50) ~ 1e-19, below the double
// rounding of a log-likelihood.
constexpr double kWindowCut = 50.0;

bool isPositiveFinite(double value) {
    return std::isfinite(value) && value > 0.0;
}

} // namespace

CaseAModel::CaseAModel(Grid En_grid, Grid theta_grid, std::vector<double> log_prior_theta,
                       CaseAResolution resolution)
    : En_grid_(std::move(En_grid)), theta_grid_(std::move(theta_grid)),
      log_prior_theta_(std::move(log_prior_theta)), resolution_(resolution),
      inv_sigma_Ep_(1.0 / resolution.sigma_Ep), inv_sigma_theta_(1.0 / resolution.sigma_theta) {
    if (!isPositiveFinite(resolution_.sigma_Ep) || !isPositiveFinite(resolution_.sigma_theta)) {
        throw std::invalid_argument("Case A resolutions must be positive and finite");
    }
    if (theta_grid_.front() < 0.0 || theta_grid_.back() > kHalf * units::pi) {
        throw std::invalid_argument("the theta_p grid must lie in [0, pi/2]");
    }
    if (log_prior_theta_.size() != theta_grid_.size()) {
        throw std::invalid_argument("the theta_p prior does not match the theta_p grid");
    }
    for (const double value : log_prior_theta_) {
        if (std::isnan(value) || (std::isinf(value) && value > 0.0)) {
            throw std::invalid_argument("the theta_p prior has a NaN or +inf term");
        }
    }
    if (!(En_grid_.front() > 0.0)) {
        throw std::invalid_argument("the E_n grid must be positive");
    }
    const std::size_t n_theta = theta_grid_.size();
    Ep_table_.resize(En_grid_.size() * n_theta);
    for (std::size_t i = 0; i < En_grid_.size(); ++i) {
        for (std::size_t j = 0; j < n_theta; ++j) {
            Ep_table_[(i * n_theta) + j] = kinematics::protonEnergy(En_grid_[i], theta_grid_[j]);
        }
    }
    // The window of rowLogSumExp needs every row non-increasing in theta_j;
    // the exact kinematics is, but it is checked rather than assumed.
    for (std::size_t i = 0; i < En_grid_.size() && Ep_decreasing_; ++i) {
        const double* row = Ep_table_.data() + (i * n_theta);
        Ep_decreasing_ = std::is_sorted(row, row + n_theta, std::greater<>());
    }
    max_log_prior_theta_ = *std::max_element(log_prior_theta_.begin(), log_prior_theta_.end());
}

double CaseAModel::rowLogSumExp(double Ep_hat, double theta_hat, const double* Ep_row,
                                std::size_t j_theta, Workspace& work) const {
    const std::size_t n_theta = theta_grid_.size();
    const auto term = [&](std::size_t j) {
        const double z = (Ep_hat - Ep_row[j]) * inv_sigma_Ep_;
        return (-kHalf * z * z) + work.angular[j];
    };
    std::size_t begin = 0;
    std::size_t end = n_theta;
    if (Ep_decreasing_) {
        // Lower bound L on the row maximum: the terms at the cells nearest to
        // Ep_hat (along the row) and to theta_hat.
        const double* Ep_end = Ep_row + n_theta;
        const auto nearest = static_cast<std::size_t>(
            std::lower_bound(Ep_row, Ep_end, Ep_hat, std::greater<>()) - Ep_row);
        double lower = term(j_theta);
        if (nearest < n_theta) {
            lower = std::max(lower, term(nearest));
        }
        if (nearest > 0) {
            lower = std::max(lower, term(nearest - 1));
        }
        if (std::isfinite(lower)) {
            // A term reaching L - C needs both Gaussian exponents within
            // a_max + C - L of 0, since the other part is <= a_max.
            const double radius = std::sqrt(2.0 * (max_log_prior_theta_ + kWindowCut - lower));
            const double Ep_hi = Ep_hat + (radius * resolution_.sigma_Ep);
            const double Ep_lo = Ep_hat - (radius * resolution_.sigma_Ep);
            const double theta_lo = theta_hat - (radius * resolution_.sigma_theta);
            const double theta_hi = theta_hat + (radius * resolution_.sigma_theta);
            const auto& theta = theta_grid_.values();
            const auto first_Ep = static_cast<std::size_t>(
                std::lower_bound(Ep_row, Ep_end, Ep_hi, std::greater<>()) - Ep_row);
            const auto last_Ep = static_cast<std::size_t>(
                std::upper_bound(Ep_row, Ep_end, Ep_lo, std::greater<>()) - Ep_row);
            const auto first_theta = static_cast<std::size_t>(
                std::lower_bound(theta.begin(), theta.end(), theta_lo) - theta.begin());
            const auto last_theta = static_cast<std::size_t>(
                std::upper_bound(theta.begin(), theta.end(), theta_hi) - theta.begin());
            begin = std::max(first_Ep, first_theta);
            end = std::min(last_Ep, last_theta);
            if (begin >= end) {
                // Unreachable in exact arithmetic (the cells giving L are in
                // the window); kept as a safe fallback.
                begin = 0;
                end = n_theta;
            }
        }
    }
    for (std::size_t j = begin; j < end; ++j) {
        work.terms[j - begin] = term(j);
    }
    // Every term is finite or -inf: the observables were checked by the
    // caller and the prior and the kinematic table at construction.
    return logSumExpUnchecked(work.terms.data(), end - begin);
}

void CaseAModel::fillLogLikelihood(double Ep_hat, double theta_hat, double* out,
                                   Workspace& work) const {
    if (!std::isfinite(Ep_hat) || !std::isfinite(theta_hat)) {
        throw std::invalid_argument("Case A observables must be finite");
    }
    const std::size_t n_theta = theta_grid_.size();
    // The angular term and the prior do not depend on E_n: once per event.
    work.angular.resize(n_theta);
    for (std::size_t j = 0; j < n_theta; ++j) {
        const double z = (theta_hat - theta_grid_[j]) * inv_sigma_theta_;
        work.angular[j] = (-kHalf * z * z) + log_prior_theta_[j];
    }
    work.terms.resize(n_theta);
    // The theta cell nearest to theta_hat, a candidate for the row maximum.
    const auto& theta = theta_grid_.values();
    std::size_t j_theta = static_cast<std::size_t>(
        std::lower_bound(theta.begin(), theta.end(), theta_hat) - theta.begin());
    if (j_theta == n_theta) {
        j_theta = n_theta - 1;
    } else if (j_theta > 0 && (theta_hat - theta[j_theta - 1]) < (theta[j_theta] - theta_hat)) {
        --j_theta;
    }
    for (std::size_t i = 0; i < En_grid_.size(); ++i) {
        out[i] = rowLogSumExp(Ep_hat, theta_hat, Ep_table_.data() + (i * n_theta), j_theta, work);
    }
}

std::vector<double> CaseAModel::logLikelihood(double Ep_hat, double theta_hat) const {
    std::vector<double> result(En_grid_.size());
    Workspace work;
    fillLogLikelihood(Ep_hat, theta_hat, result.data(), work);
    return result;
}

LogLikelihoodMatrix CaseAModel::logLikelihoods(const std::vector<double>& Ep_hat,
                                               const std::vector<double>& theta_hat) const {
    if (Ep_hat.size() != theta_hat.size()) {
        throw std::invalid_argument("Ep_hat and theta_hat have different lengths");
    }
    LogLikelihoodMatrix matrix;
    matrix.events = Ep_hat.size();
    matrix.points = En_grid_.size();
    matrix.values.resize(matrix.events * matrix.points);
    Workspace work;
    for (std::size_t k = 0; k < matrix.events; ++k) {
        fillLogLikelihood(Ep_hat[k], theta_hat[k], matrix.values.data() + (k * matrix.points),
                          work);
    }
    return matrix;
}

LogLikelihoodMatrix CaseAModel::logLikelihoods(const EventSample& sample) const {
    return logLikelihoods(sample.reco_Ep, sample.reco_theta_p);
}

std::vector<double> singleEventLogPosterior(const CaseAModel& model, double Ep_hat,
                                            double theta_hat,
                                            const std::vector<double>& log_prior_En) {
    if (log_prior_En.size() != model.energyGrid().size()) {
        throw std::invalid_argument("the E_n prior does not match the E_n grid");
    }
    std::vector<double> result = model.logLikelihood(Ep_hat, theta_hat);
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] += log_prior_En[i];
    }
    return result;
}

} // namespace riptide::stats
