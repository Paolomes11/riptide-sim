#pragma once

// Case A (D33): one shared neutron energy E_n, known direction (the beam
// axis). Each event gives (reco_Ep, reco_theta_p); its likelihood on the E_n
// grid marginalises the true recoil angle theta_p over a grid on [0, pi/2]:
//
//   loglik(E_n,i) = log sum_j exp[ -1/2 ((Ep_hat - Ep(E_n,i, theta_j)) / sigma_Ep)^2
//                                  -1/2 ((theta_hat - theta_j) / sigma_theta)^2
//                                  + log prior(theta_j) ]
//
// with Ep(E_n, theta_p) the exact n-p kinematics (D32). The Gaussian
// normalisations do not depend on E_n and are left out, so the likelihood is
// defined up to one constant per event.

#include <cstddef>
#include <vector>

#include "riptide/common/event_sample.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

struct CaseAResolution {
    double sigma_Ep = 0.0;    // MeV
    double sigma_theta = 0.0; // rad
};

// Per-event log-likelihoods, row-major: row k is event k on the E_n grid.
struct LogLikelihoodMatrix {
    std::size_t events = 0;
    std::size_t points = 0;
    std::vector<double> values;

    [[nodiscard]] const double* row(std::size_t k) const noexcept {
        return values.data() + (k * points);
    }
};

class CaseAModel {
public:
    // log_prior_theta holds the log mass of each theta cell (priors.hpp).
    // Throws std::invalid_argument for non-positive or non-finite
    // resolutions, a theta grid outside [0, pi/2], a prior of the wrong size
    // or a prior term that is NaN or +inf.
    CaseAModel(Grid En_grid, Grid theta_grid, std::vector<double> log_prior_theta,
               CaseAResolution resolution);

    [[nodiscard]] const Grid& energyGrid() const noexcept {
        return En_grid_;
    }
    [[nodiscard]] const Grid& recoilAngleGrid() const noexcept {
        return theta_grid_;
    }
    [[nodiscard]] const CaseAResolution& resolution() const noexcept {
        return resolution_;
    }

    // One event. Ep_hat may be negative (the smearing of a soft proton);
    // both observables must be finite, else std::invalid_argument.
    [[nodiscard]] std::vector<double> logLikelihood(double Ep_hat, double theta_hat) const;

    [[nodiscard]] LogLikelihoodMatrix logLikelihoods(const std::vector<double>& Ep_hat,
                                                     const std::vector<double>& theta_hat) const;

    // reco_Ep and reco_theta_p of a prepared sample.
    [[nodiscard]] LogLikelihoodMatrix logLikelihoods(const EventSample& sample) const;

private:
    // Buffers reused across the events of one call, so the hot loop does not
    // allocate.
    struct Workspace {
        std::vector<double> angular; // theta term + log prior, per theta_j
        std::vector<double> terms;   // the summands of one E_n row
    };

    void fillLogLikelihood(double Ep_hat, double theta_hat, double* out, Workspace& work) const;

    // logsumexp of one E_n row over the theta_p cells that can matter. A
    // lower bound L on the row maximum bounds every term of a cell outside a
    // window in Ep and in theta below L - C; those cells are skipped.
    [[nodiscard]] double rowLogSumExp(double Ep_hat, double theta_hat, const double* Ep_row,
                                      std::size_t j_theta, Workspace& work) const;

    Grid En_grid_;
    Grid theta_grid_;
    std::vector<double> log_prior_theta_;
    CaseAResolution resolution_;
    // Multiplying by the inverse is cheaper than dividing in the hot loop; the
    // difference is at the last bit of z.
    double inv_sigma_Ep_;
    double inv_sigma_theta_;
    // Ep(E_n,i, theta_j) at [i * n_theta + j]: the kinematics is evaluated once
    // per model, not once per event.
    std::vector<double> Ep_table_;
    // Every row of Ep_table_ is non-increasing in theta_j, so the Ep window is
    // a contiguous range found by bisection; otherwise rows are summed whole.
    bool Ep_decreasing_ = true;
    double max_log_prior_theta_ = 0.0;
};

// loglik + log prior on the E_n grid: the unnormalised posterior of one event.
[[nodiscard]] std::vector<double> singleEventLogPosterior(const CaseAModel& model, double Ep_hat,
                                                          double theta_hat,
                                                          const std::vector<double>& log_prior_En);

} // namespace riptide::stats
