#include "riptide/stats/case_b.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/log_sum_exp.hpp"
#include "riptide/stats/vmf_kernel.hpp"

namespace riptide::stats {

namespace {

constexpr double kHalf = 0.5;

// Cut C of the E_n window, as in Case A (D35 addendum): every skipped term is
// at least C below the row maximum, so the skipped mass is at most
// n_En * exp(-C) relative to the kept one: 500 * exp(-50) ~ 1e-19, below the
// double rounding of a log-likelihood.
constexpr double kWindowCut = 50.0;

// The interpolation position stays this far below the last cell, so that the
// upper neighbour exists; the same clamp as the toy model.
constexpr double kLastCellMargin = 1e-9;

bool isPositiveFinite(double value) {
    return std::isfinite(value) && value > 0.0;
}

} // namespace

double logRecoilDensityRatio(double theta_p, double En) noexcept {
    const double two_m = 2.0 * kinematics::kProtonMass;
    const double sin_theta = std::sin(theta_p);
    return std::log(two_m * (two_m + En)) - (2.0 * std::log(two_m + (En * sin_theta * sin_theta)));
}

std::vector<double> logTrackKernelTable(const Grid& theta_obs_grid, const Grid& theta_grid,
                                        double sigma_theta) {
    const double kappa = vmfConcentration(sigma_theta);
    const std::size_t n_theta = theta_grid.size();
    // Trapezoid weights: the end points count half.
    std::vector<double> log_weight(n_theta);
    for (std::size_t j = 0; j < n_theta; ++j) {
        const double weight = (j == 0 || j + 1 == n_theta) ? kHalf : 1.0;
        log_weight[j] = std::log(weight * theta_grid.step()) + logTrackDensity(theta_grid[j]);
    }
    std::vector<double> table(theta_obs_grid.size() * n_theta);
    for (std::size_t o = 0; o < theta_obs_grid.size(); ++o) {
        for (std::size_t j = 0; j < n_theta; ++j) {
            table[(o * n_theta) + j] =
                logMeridianVmfKernel(theta_obs_grid[o], theta_grid[j], kappa) + log_weight[j];
        }
    }
    return table;
}

CaseBModel::CaseBModel(Grid En_grid, Grid theta_grid, Grid theta_obs_grid,
                       std::vector<double> log_prior_En, CaseBResolution resolution,
                       const std::vector<double>& log_acceptance)
    : En_grid_(std::move(En_grid)), theta_grid_(std::move(theta_grid)),
      theta_obs_grid_(std::move(theta_obs_grid)), log_prior_En_(std::move(log_prior_En)),
      resolution_(resolution), inv_sigma_Ep_(1.0 / resolution.sigma_Ep),
      inv_obs_step_(1.0 / theta_obs_grid_.step()) {
    if (!isPositiveFinite(resolution_.sigma_Ep) || !isPositiveFinite(resolution_.sigma_theta)) {
        throw std::invalid_argument("Case B resolutions must be positive and finite");
    }
    if (theta_grid_.front() < 0.0 || theta_grid_.back() > kHalf * units::pi) {
        throw std::invalid_argument("the theta_p grid must lie in [0, pi/2]");
    }
    // The interpolation maps theta_obs to a position by the step alone, so
    // the grid must start at 0 and end at pi.
    if (std::abs(theta_obs_grid_.front()) > 0.0 ||
        std::abs(theta_obs_grid_.back() - units::pi) > 0.0) {
        throw std::invalid_argument("the theta_obs grid must be [0, pi]");
    }
    if (!(En_grid_.front() > 0.0)) {
        throw std::invalid_argument("the E_n grid must be positive");
    }
    const std::size_t n_En = En_grid_.size();
    if (log_prior_En_.size() != n_En) {
        throw std::invalid_argument("the E_n prior does not match the E_n grid");
    }
    bool any_finite = false;
    for (const double value : log_prior_En_) {
        if (std::isnan(value) || (std::isinf(value) && value > 0.0)) {
            throw std::invalid_argument("the E_n prior has a NaN or +inf term");
        }
        any_finite = any_finite || std::isfinite(value);
    }
    if (!any_finite) {
        throw std::invalid_argument("the E_n prior has no finite term");
    }
    const std::size_t n_theta = theta_grid_.size();
    if (!log_acceptance.empty()) {
        if (log_acceptance.size() != n_theta * n_En) {
            throw std::invalid_argument("the acceptance table does not match the grids");
        }
        if (!std::all_of(log_acceptance.begin(), log_acceptance.end(),
                         [](double value) { return std::isfinite(value); })) {
            throw std::invalid_argument("the acceptance table has a non-finite term");
        }
        log_acceptance_ = log_acceptance;
    }
    Ep_table_.resize(n_theta * n_En);
    energy_weight_.resize(n_theta * n_En);
    max_energy_weight_.resize(n_theta);
    for (std::size_t j = 0; j < n_theta; ++j) {
        double row_max = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < n_En; ++i) {
            const std::size_t index = (j * n_En) + i;
            Ep_table_[index] = kinematics::protonEnergy(En_grid_[i], theta_grid_[j]);
            energy_weight_[index] =
                logRecoilDensityRatio(theta_grid_[j], En_grid_[i]) + log_prior_En_[i];
            if (!log_acceptance_.empty()) {
                energy_weight_[index] += log_acceptance_[index];
            }
            row_max = std::max(row_max, energy_weight_[index]);
        }
        max_energy_weight_[j] = row_max;
        // The window needs every row non-decreasing in E_n; the exact
        // kinematics is, but it is checked rather than assumed.
        const double* row = Ep_table_.data() + (j * n_En);
        Ep_increasing_ = Ep_increasing_ && std::is_sorted(row, row + n_En);
    }
    kernel_ = logTrackKernelTable(theta_obs_grid_, theta_grid_, resolution_.sigma_theta);
}

double CaseBModel::energyLogSumExp(double Ep_hat, std::size_t j, Workspace& work) const {
    const std::size_t n_En = En_grid_.size();
    const double* Ep_row = Ep_table_.data() + (j * n_En);
    const double* weight_row = energy_weight_.data() + (j * n_En);
    const auto term = [&](std::size_t i) {
        const double z = (Ep_hat - Ep_row[i]) * inv_sigma_Ep_;
        return (-kHalf * z * z) + weight_row[i];
    };
    std::size_t begin = 0;
    std::size_t end = n_En;
    if (Ep_increasing_) {
        const double* Ep_end = Ep_row + n_En;
        const auto nearest =
            static_cast<std::size_t>(std::lower_bound(Ep_row, Ep_end, Ep_hat) - Ep_row);
        double lower = -std::numeric_limits<double>::infinity();
        if (nearest < n_En) {
            lower = term(nearest);
        }
        if (nearest > 0) {
            lower = std::max(lower, term(nearest - 1));
        }
        if (std::isfinite(lower)) {
            // A term reaching L - C needs its Gaussian exponent within
            // w_max + C - L of 0, since its weight is <= w_max.
            const double radius = std::sqrt(2.0 * (max_energy_weight_[j] + kWindowCut - lower));
            const double Ep_lo = Ep_hat - (radius * resolution_.sigma_Ep);
            const double Ep_hi = Ep_hat + (radius * resolution_.sigma_Ep);
            begin = static_cast<std::size_t>(std::lower_bound(Ep_row, Ep_end, Ep_lo) - Ep_row);
            end = static_cast<std::size_t>(std::upper_bound(Ep_row, Ep_end, Ep_hi) - Ep_row);
            if (begin >= end) {
                // Unreachable in exact arithmetic (the cells giving L are in
                // the window); kept as a safe fallback.
                begin = 0;
                end = n_En;
            }
        }
    }
    for (std::size_t i = begin; i < end; ++i) {
        work.terms[i - begin] = term(i);
    }
    // Every term is finite or -inf: Ep_hat was checked by the caller, the
    // prior and the tables at construction.
    return logSumExpUnchecked(work.terms.data(), end - begin);
}

void CaseBModel::fillTrackEnergyRow(double Ep_hat, double* out, Workspace& work) const {
    if (!std::isfinite(Ep_hat)) {
        throw std::invalid_argument("Case B observables must be finite");
    }
    const std::size_t n_theta = theta_grid_.size();
    work.energy.resize(n_theta);
    work.terms.resize(std::max(n_theta, En_grid_.size()));
    for (std::size_t j = 0; j < n_theta; ++j) {
        work.energy[j] = energyLogSumExp(Ep_hat, j, work);
    }
    for (std::size_t o = 0; o < theta_obs_grid_.size(); ++o) {
        const double* kernel_row = kernel_.data() + (o * n_theta);
        for (std::size_t j = 0; j < n_theta; ++j) {
            work.terms[j] = kernel_row[j] + work.energy[j];
        }
        out[o] = logSumExpUnchecked(work.terms.data(), n_theta);
    }
}

std::vector<double> CaseBModel::trackEnergyRow(double Ep_hat) const {
    std::vector<double> result(theta_obs_grid_.size());
    Workspace work;
    fillTrackEnergyRow(Ep_hat, result.data(), work);
    return result;
}

TrackEnergyTable CaseBModel::trackEnergyTable(const std::vector<double>& Ep_hat) const {
    TrackEnergyTable table;
    table.events = Ep_hat.size();
    table.points = theta_obs_grid_.size();
    table.values.resize(table.events * table.points);
    Workspace work;
    for (std::size_t k = 0; k < table.events; ++k) {
        fillTrackEnergyRow(Ep_hat[k], table.values.data() + (k * table.points), work);
    }
    return table;
}

double CaseBModel::interpolate(const double* row, double theta_obs) const noexcept {
    const auto last = static_cast<double>(theta_obs_grid_.size() - 1);
    const double position = std::clamp(theta_obs * inv_obs_step_, 0.0, last - kLastCellMargin);
    const auto index = static_cast<std::size_t>(position);
    const double fraction = position - static_cast<double>(index);
    return row[index] + (fraction * (row[index + 1] - row[index]));
}

std::vector<double>
CaseBModel::summedLogLikelihood(const TrackEnergyTable& table,
                                const std::vector<geometry::Direction>& tracks,
                                const std::vector<geometry::Direction>& candidates) const {
    if (tracks.size() != table.events) {
        throw std::invalid_argument("the tracks do not match the rows of the table");
    }
    if (table.points != theta_obs_grid_.size()) {
        throw std::invalid_argument("the table does not match the theta_obs grid");
    }
    std::vector<double> result(candidates.size(), 0.0);
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        // Events summed in a fixed order: the result is the same bit for bit
        // whatever splits the candidates.
        double sum = 0.0;
        for (std::size_t k = 0; k < tracks.size(); ++k) {
            sum += interpolate(table.row(k), geometry::angleBetween(tracks[k], candidates[c]));
        }
        result[c] = sum;
    }
    return result;
}

std::vector<double> CaseBModel::energyLogPosterior(double Ep_hat, double theta_obs) const {
    if (!std::isfinite(Ep_hat) || !std::isfinite(theta_obs)) {
        throw std::invalid_argument("Case B observables must be finite");
    }
    const std::size_t n_theta = theta_grid_.size();
    const std::size_t n_En = En_grid_.size();
    const double kappa = vmfConcentration(resolution_.sigma_theta);
    std::vector<double> angular(n_theta);
    for (std::size_t j = 0; j < n_theta; ++j) {
        const double weight = (j == 0 || j + 1 == n_theta) ? kHalf : 1.0;
        angular[j] = logMeridianVmfKernel(theta_obs, theta_grid_[j], kappa) +
                     logTrackDensity(theta_grid_[j]) + std::log(weight * theta_grid_.step());
    }
    std::vector<double> result(n_En);
    std::vector<double> terms(n_theta);
    for (std::size_t i = 0; i < n_En; ++i) {
        for (std::size_t j = 0; j < n_theta; ++j) {
            const std::size_t index = (j * n_En) + i;
            const double z = (Ep_hat - Ep_table_[index]) * inv_sigma_Ep_;
            const double acceptance = log_acceptance_.empty() ? 0.0 : log_acceptance_[index];
            terms[j] = angular[j] + (-kHalf * z * z) +
                       logRecoilDensityRatio(theta_grid_[j], En_grid_[i]) + acceptance;
        }
        result[i] = logSumExpUnchecked(terms.data(), n_theta) + log_prior_En_[i];
    }
    return result;
}

} // namespace riptide::stats
