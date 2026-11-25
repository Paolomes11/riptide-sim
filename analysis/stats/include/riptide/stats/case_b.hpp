#pragma once

// Case B (D41): one shared, unknown neutron direction Omega_n; E_n is a
// nuisance of each event with a proper prior. Each event gives (reco_Ep,
// track); its likelihood for a candidate Omega depends on Omega only through
// the angle theta_obs between the observed track and Omega:
//
//   g_k(theta_j) = log sum_i exp[ -1/2 ((Ep_hat - Ep(E_n,i, theta_j)) / sigma_Ep)^2
//                                 + log rho(theta_j, E_n,i) + log a(theta_j, E_n,i)
//                                 + log prior(E_n,i) ]
//   h_k(theta_obs) = log sum_j exp[ log K(theta_obs | theta_j) + log(cos theta_j / pi)
//                                   + log w_j + g_k(theta_j) ]
//
// with K the azimuth-integrated vMF kernel (vmf_kernel.hpp), w_j the
// trapezoid weights of the theta_p grid and rho the relativistic factor of
// the track density and a the optional selection acceptance (1 for the ideal
// detector). h_k is tabulated once per event on a uniform theta_obs
// grid over [0, pi] and interpolated linearly at each candidate, so the cost
// per candidate does not depend on the E_n and theta_p grids. The Gaussian
// normalisations are left out: the likelihood is defined up to one constant
// per event.

#include <cstddef>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// Defaults of the toy model, so that its numbers compare at the same
// discretisation (the E_n grid is kEnergyGridPoints).
inline constexpr std::size_t kTrackAngleGridPoints = 400;
inline constexpr std::size_t kTrackObservedGridPoints = 801;

struct CaseBResolution {
    double sigma_Ep = 0.0;    // MeV
    double sigma_theta = 0.0; // rad, per component of the track offset
};

// log of the relativistic density of the recoil direction about Omega_n over
// its non-relativistic limit cos(theta_p) / pi, for equal nucleon masses m:
// log[2m (2m + E_n) / (2m + E_n sin^2 theta_p)^2]. It is 0 at E_n -> 0 and
// within ~0.011 of 0 at 20 MeV.
[[nodiscard]] double logRecoilDensityRatio(double theta_p, double En) noexcept;

// log weights of the theta_p quadrature for an observed angle, row-major
// n_obs x n_theta: log K(theta_obs | theta_j) + log(cos theta_j / pi) + log w_j.
// exp of a row summed over j is the density per steradian of the observed
// track in the non-relativistic limit. Throws std::invalid_argument for a
// non-positive or non-finite sigma_theta.
[[nodiscard]] std::vector<double> logTrackKernelTable(const Grid& theta_obs_grid,
                                                      const Grid& theta_grid, double sigma_theta);

// h_k on the theta_obs grid, row-major: row k is event k.
struct TrackEnergyTable {
    std::size_t events = 0;
    std::size_t points = 0;
    std::vector<double> values;

    [[nodiscard]] const double* row(std::size_t k) const noexcept {
        return values.data() + (k * points);
    }
};

class CaseBModel {
public:
    // theta_grid on [0, pi/2], theta_obs_grid exactly [0, pi], log_prior_En a
    // log density on En_grid (priors.hpp). log_acceptance is empty for the
    // ideal detector, or the caseBAcceptanceLogWeights table (acceptance.hpp)
    // on the same grids, added to the recoil density. Throws
    // std::invalid_argument for non-positive or non-finite resolutions, grids
    // outside those ranges, an E_n grid that is not positive, a prior of the
    // wrong size, with a NaN or +inf term or with no finite term, and an
    // acceptance table of the wrong size or with a non-finite term.
    CaseBModel(Grid En_grid, Grid theta_grid, Grid theta_obs_grid, std::vector<double> log_prior_En,
               CaseBResolution resolution, const std::vector<double>& log_acceptance = {});

    [[nodiscard]] const Grid& energyGrid() const noexcept {
        return En_grid_;
    }
    [[nodiscard]] const Grid& recoilAngleGrid() const noexcept {
        return theta_grid_;
    }
    [[nodiscard]] const Grid& observedAngleGrid() const noexcept {
        return theta_obs_grid_;
    }
    [[nodiscard]] const CaseBResolution& resolution() const noexcept {
        return resolution_;
    }
    [[nodiscard]] bool hasAcceptance() const noexcept {
        return !log_acceptance_.empty();
    }

    // h_k of one event on the theta_obs grid. Ep_hat must be finite (it may
    // be negative), else std::invalid_argument.
    [[nodiscard]] std::vector<double> trackEnergyRow(double Ep_hat) const;

    [[nodiscard]] TrackEnergyTable trackEnergyTable(const std::vector<double>& Ep_hat) const;

    // Linear interpolation of a row of the table at theta_obs in [0, pi]; the
    // position is clamped just below the last cell, as in the toy model.
    [[nodiscard]] double interpolate(const double* row, double theta_obs) const noexcept;

    // sum_k h_k(angle(track_k, candidate_c)) for each candidate: the log-
    // likelihood of the sample on a direction grid, prior not included.
    // tracks and the rows of the table must match in number, else
    // std::invalid_argument. Unit vectors are expected.
    [[nodiscard]] std::vector<double>
    summedLogLikelihood(const TrackEnergyTable& table,
                        const std::vector<geometry::Direction>& tracks,
                        const std::vector<geometry::Direction>& candidates) const;

    // Unnormalised log posterior of E_n for one event at a fixed observed
    // angle: log prior(E_n,i) + log sum_j exp[kernel + rho + Gaussian]. The
    // kernel is evaluated at theta_obs itself, not on the table, and no window
    // is used: a diagnostic, and the Omega_n || z limit against Case A.
    [[nodiscard]] std::vector<double> energyLogPosterior(double Ep_hat, double theta_obs) const;

private:
    // Buffers reused across the events of one call, so the loop does not
    // allocate.
    struct Workspace {
        std::vector<double> energy; // g_k(theta_j)
        std::vector<double> terms;  // the summands of one logsumexp
    };

    void fillTrackEnergyRow(double Ep_hat, double* out, Workspace& work) const;

    // g_k(theta_j) over the E_n cells that can matter: a lower bound L on the
    // row maximum bounds every term outside a window in Ep below L - C.
    [[nodiscard]] double energyLogSumExp(double Ep_hat, std::size_t j, Workspace& work) const;

    Grid En_grid_;
    Grid theta_grid_;
    Grid theta_obs_grid_;
    std::vector<double> log_prior_En_;
    CaseBResolution resolution_;
    double inv_sigma_Ep_;
    double inv_obs_step_;
    // Ep(E_n,i, theta_j) at [j * n_En + i]: along a row Ep grows with E_n.
    std::vector<double> Ep_table_;
    // log rho(theta_j, E_n,i) + log prior(E_n,i), plus the acceptance term
    // when given, at the same index.
    std::vector<double> energy_weight_;
    // The acceptance table as given, empty for the ideal detector.
    std::vector<double> log_acceptance_;
    // Maximum of energy_weight_ per theta_j row, for the window bound.
    std::vector<double> max_energy_weight_;
    // Every row of Ep_table_ is non-decreasing in E_n, so the window is found
    // by bisection; otherwise rows are summed whole.
    bool Ep_increasing_ = true;
    // logTrackKernelTable, n_obs x n_theta.
    std::vector<double> kernel_;
};

} // namespace riptide::stats
