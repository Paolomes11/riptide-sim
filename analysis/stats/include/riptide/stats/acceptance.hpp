#pragma once

// theta_p acceptance of the Case A selection (D40). The selected events of a
// simulated run, histogrammed in (E_n, theta_p), give the distribution of the
// true recoil angle at each E_n: it folds together the n-p angular
// distribution and the selection (single n-p scatter, proton contained).
// Dividing by the isotropic distribution leaves a smooth ratio r(E_n, theta_p);
// the model prior is then q(theta | E_n) proportional to r * pi_iso, with
// pi_iso the exact relativistic isotropic distribution, so the isotropic case
// r = 1 is reproduced exactly.
//
// Only the shape in theta_p at fixed E_n enters: Case A conditions on the
// selection, so the overall acceptance A(E_n) cancels in every row.

#include <cstddef>
#include <vector>

#include "riptide/stats/grid.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// Counts of selected events. counts[k * n_theta + t] is E_n bin k
// [En_edges[k], En_edges[k+1]) and theta_p bin t [theta_edges[t],
// theta_edges[t+1]), with n_theta = theta_edges.size() - 1. Edges in MeV and
// rad.
struct AcceptanceHistogram {
    std::vector<double> En_edges;
    std::vector<double> theta_edges;
    std::vector<double> counts;
};

// Probability that an isotropic (in the CM) n-p scatter at E_n gives a proton
// at a lab angle in [lo, hi], 0 <= lo <= hi <= pi/2. The recoil energy is
// linear in the CM cosine for any masses, so Ep is uniform on [0, Ep_max] and
// the mass is u(lo) - u(hi) with u = Ep(E_n, theta) / Ep_max(E_n), exact and
// relativistic. Throws std::invalid_argument for a bad interval or E_n <= 0.
[[nodiscard]] double isotropicRecoilAngleMass(double En, double lo, double hi);

// r(E_n, theta_p) = (n_kt / N_k) / m_iso(t; E_k): the selected fraction of
// bin t in E_n bin k over the isotropic mass of that bin, at the bin centre
// E_k. Evaluated by bilinear interpolation between bin centres and held
// constant beyond the first and last centre in each variable.
class AcceptanceRatio {
public:
    // Throws std::invalid_argument unless both edge lists are finite and
    // strictly increasing with at least one bin, E_n edges are positive,
    // theta edges lie in [0, pi/2], the counts have one entry per bin and
    // every count is finite and positive (an empty bin would put a zero in
    // the prior and is a sign of too fine a binning).
    explicit AcceptanceRatio(const AcceptanceHistogram& histogram);

    [[nodiscard]] double operator()(double En, double theta) const noexcept;

    // The largest value of r anywhere: bilinear interpolation never exceeds
    // the largest node, so it is the maximum over the bins.
    [[nodiscard]] double maximum() const noexcept {
        return maximum_;
    }

    [[nodiscard]] const std::vector<double>& energyCentres() const noexcept {
        return En_centres_;
    }
    [[nodiscard]] const std::vector<double>& angleCentres() const noexcept {
        return theta_centres_;
    }

private:
    std::vector<double> En_centres_;
    std::vector<double> theta_centres_;
    std::vector<double> ratio_; // [k * n_theta + t]
    double maximum_ = 0.0;
};

// Log prior of theta_p per E_n row, for the per-row CaseAModel constructor:
// row i, cell j is log of the normalised integral of r(E_n,i, theta) over the
// isotropic distribution on the theta cell j ([theta_j - step/2,
// theta_j + step/2] clipped to [0, pi/2], as physicalRecoilAngleLogPrior).
// Each row sums to 1 up to rounding.
[[nodiscard]] std::vector<std::vector<double>>
acceptanceRecoilAngleLogPrior(const AcceptanceRatio& ratio, const Grid& En_grid,
                              const Grid& theta_grid);

} // namespace riptide::stats
