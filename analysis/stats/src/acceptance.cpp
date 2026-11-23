#include "riptide/stats/acceptance.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"

namespace riptide::stats {

namespace {

constexpr double kHalf = 0.5;
constexpr double kThetaMax = kHalf * units::pi;

// Sub-intervals per theta cell, equal in u = Ep / Ep_max, of the midpoint
// rule for the cell mass of r * pi_iso. pi_iso is uniform in u, so the rule
// is exact for r constant on the cell; r varies on the scale of a 5 deg bin,
// several model cells, so a few sub-intervals suffice.
constexpr std::size_t kCellSubdivisions = 8;

bool strictlyIncreasingFinite(const std::vector<double>& edges) {
    if (edges.size() < 2) {
        return false;
    }
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (!std::isfinite(edges[i]) || (i > 0 && !(edges[i] > edges[i - 1]))) {
            return false;
        }
    }
    return true;
}

std::vector<double> centres(const std::vector<double>& edges) {
    std::vector<double> result(edges.size() - 1);
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = kHalf * (edges[i] + edges[i + 1]);
    }
    return result;
}

// Lower node and weight of the upper node for x among increasing centres,
// with x clamped to [front, back]: constant extrapolation. One centre gives
// weight 0 on node 0.
struct Bracket {
    std::size_t lower = 0;
    double weight = 0.0;
};

Bracket bracket(const std::vector<double>& nodes, double x) {
    if (nodes.size() == 1 || !(x > nodes.front())) {
        return {0, 0.0};
    }
    if (!(x < nodes.back())) {
        return {nodes.size() - 2, 1.0};
    }
    const auto upper =
        static_cast<std::size_t>(std::upper_bound(nodes.begin(), nodes.end(), x) - nodes.begin());
    const std::size_t lower = upper - 1;
    return {lower, (x - nodes[lower]) / (nodes[upper] - nodes[lower])};
}

// u(theta) = Ep(E_n, theta) / Ep_max(E_n), 1 at theta = 0 and 0 at pi/2.
double energyFraction(double En, double theta) {
    return kinematics::protonEnergy(En, theta) / kinematics::maxProtonEnergy(En);
}

} // namespace

double isotropicRecoilAngleMass(double En, double lo, double hi) {
    if (!(En > 0.0) || !std::isfinite(En) || !(lo >= 0.0) || !(hi >= lo) || !(hi <= kThetaMax)) {
        throw std::invalid_argument(
            "isotropic theta_p mass needs E_n > 0 and 0 <= lo <= hi <= pi/2");
    }
    return energyFraction(En, lo) - energyFraction(En, hi);
}

AcceptanceRatio::AcceptanceRatio(const AcceptanceHistogram& histogram) {
    const auto& En_edges = histogram.En_edges;
    const auto& theta_edges = histogram.theta_edges;
    if (!strictlyIncreasingFinite(En_edges) || !(En_edges.front() > 0.0)) {
        throw std::invalid_argument("acceptance E_n edges must be positive and increasing");
    }
    if (!strictlyIncreasingFinite(theta_edges) || theta_edges.front() < 0.0 ||
        theta_edges.back() > kThetaMax) {
        throw std::invalid_argument("acceptance theta_p edges must be increasing in [0, pi/2]");
    }
    const std::size_t n_En = En_edges.size() - 1;
    const std::size_t n_theta = theta_edges.size() - 1;
    if (histogram.counts.size() != n_En * n_theta) {
        throw std::invalid_argument("acceptance counts do not match the edges");
    }
    for (const double count : histogram.counts) {
        if (!std::isfinite(count) || !(count > 0.0)) {
            throw std::invalid_argument("acceptance counts must be finite and positive");
        }
    }
    En_centres_ = centres(En_edges);
    theta_centres_ = centres(theta_edges);
    ratio_.resize(n_En * n_theta);
    for (std::size_t k = 0; k < n_En; ++k) {
        const double* row = histogram.counts.data() + (k * n_theta);
        double total = 0.0;
        for (std::size_t t = 0; t < n_theta; ++t) {
            total += row[t];
        }
        for (std::size_t t = 0; t < n_theta; ++t) {
            const double iso =
                isotropicRecoilAngleMass(En_centres_[k], theta_edges[t], theta_edges[t + 1]);
            ratio_[(k * n_theta) + t] = (row[t] / total) / iso;
        }
    }
    maximum_ = *std::max_element(ratio_.begin(), ratio_.end());
}

double AcceptanceRatio::operator()(double En, double theta) const noexcept {
    const std::size_t n_theta = theta_centres_.size();
    const Bracket e = bracket(En_centres_, En);
    const Bracket a = bracket(theta_centres_, theta);
    const std::size_t e_hi = std::min(e.lower + 1, En_centres_.size() - 1);
    const std::size_t a_hi = std::min(a.lower + 1, n_theta - 1);
    const auto at = [&](std::size_t k, std::size_t t) { return ratio_[(k * n_theta) + t]; };
    const double low_row =
        ((1.0 - a.weight) * at(e.lower, a.lower)) + (a.weight * at(e.lower, a_hi));
    const double high_row = ((1.0 - a.weight) * at(e_hi, a.lower)) + (a.weight * at(e_hi, a_hi));
    return ((1.0 - e.weight) * low_row) + (e.weight * high_row);
}

std::vector<std::vector<double>> acceptanceRecoilAngleLogPrior(const AcceptanceRatio& ratio,
                                                               const Grid& En_grid,
                                                               const Grid& theta_grid) {
    if (!(En_grid.front() > 0.0)) {
        throw std::invalid_argument("the E_n grid must be positive");
    }
    if (theta_grid.front() < 0.0 || theta_grid.back() > kThetaMax) {
        throw std::invalid_argument("the theta_p grid must lie in [0, pi/2]");
    }
    const double half = kHalf * theta_grid.step();
    std::vector<std::vector<double>> rows(En_grid.size(), std::vector<double>(theta_grid.size()));
    std::vector<double> mass(theta_grid.size());
    for (std::size_t i = 0; i < En_grid.size(); ++i) {
        const double En = En_grid[i];
        const double Ep_max = kinematics::maxProtonEnergy(En);
        double total = 0.0;
        for (std::size_t j = 0; j < theta_grid.size(); ++j) {
            const double lo = std::clamp(theta_grid[j] - half, 0.0, kThetaMax);
            const double hi = std::clamp(theta_grid[j] + half, 0.0, kThetaMax);
            // u decreases with theta: the cell is [u(hi), u(lo)] in u.
            const double u_lo = energyFraction(En, hi);
            const double du =
                (energyFraction(En, lo) - u_lo) / static_cast<double>(kCellSubdivisions);
            double sum = 0.0;
            for (std::size_t s = 0; s < kCellSubdivisions; ++s) {
                const double u_mid = u_lo + ((static_cast<double>(s) + kHalf) * du);
                sum += ratio(En, kinematics::protonAngleFromEnergy(u_mid * Ep_max, En));
            }
            mass[j] = sum * du;
            total += mass[j];
        }
        for (std::size_t j = 0; j < theta_grid.size(); ++j) {
            rows[i][j] = std::log(mass[j] / total);
        }
    }
    return rows;
}

} // namespace riptide::stats
