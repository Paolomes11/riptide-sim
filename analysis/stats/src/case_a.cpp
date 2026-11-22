#include "riptide/stats/case_a.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/log_sum_exp.hpp"

namespace riptide::stats {

namespace {

constexpr double kHalf = 0.5;

bool isPositiveFinite(double value) {
    return std::isfinite(value) && value > 0.0;
}

} // namespace

CaseAModel::CaseAModel(Grid En_grid, Grid theta_grid, std::vector<double> log_prior_theta,
                       CaseAResolution resolution)
    : En_grid_(std::move(En_grid)), theta_grid_(std::move(theta_grid)),
      log_prior_theta_(std::move(log_prior_theta)), resolution_(resolution) {
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
}

void CaseAModel::fillLogLikelihood(double Ep_hat, double theta_hat, double* out,
                                   std::vector<double>& scratch) const {
    if (!std::isfinite(Ep_hat) || !std::isfinite(theta_hat)) {
        throw std::invalid_argument("Case A observables must be finite");
    }
    const std::size_t n_theta = theta_grid_.size();
    // The angular term and the prior do not depend on E_n: once per event.
    std::vector<double> angular(n_theta);
    for (std::size_t j = 0; j < n_theta; ++j) {
        const double z = (theta_hat - theta_grid_[j]) / resolution_.sigma_theta;
        angular[j] = (-kHalf * z * z) + log_prior_theta_[j];
    }
    scratch.resize(n_theta);
    for (std::size_t i = 0; i < En_grid_.size(); ++i) {
        const double* Ep_row = Ep_table_.data() + (i * n_theta);
        for (std::size_t j = 0; j < n_theta; ++j) {
            const double z = (Ep_hat - Ep_row[j]) / resolution_.sigma_Ep;
            scratch[j] = (-kHalf * z * z) + angular[j];
        }
        out[i] = logSumExp(scratch);
    }
}

std::vector<double> CaseAModel::logLikelihood(double Ep_hat, double theta_hat) const {
    std::vector<double> result(En_grid_.size());
    std::vector<double> scratch;
    fillLogLikelihood(Ep_hat, theta_hat, result.data(), scratch);
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
    std::vector<double> scratch;
    for (std::size_t k = 0; k < matrix.events; ++k) {
        fillLogLikelihood(Ep_hat[k], theta_hat[k], matrix.values.data() + (k * matrix.points),
                          scratch);
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
