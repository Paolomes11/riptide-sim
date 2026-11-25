#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/acceptance.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/log_sum_exp.hpp"
#include "riptide/stats/priors.hpp"

namespace st = riptide::stats;
namespace kin = riptide::kinematics;
namespace geo = riptide::geometry;
using Catch::Matchers::WithinAbs;
using riptide::units::pi;

namespace {

// Resolutions of the toy calibration (D43).
constexpr st::CaseBResolution kToyResolution{0.10, 0.08};
constexpr double kEnMin = 10.0;
constexpr double kEnMax = 20.0;

st::Grid observedGrid(std::size_t n = st::kTrackObservedGridPoints) {
    return {0.0, pi, n};
}

st::CaseBModel toyModel(st::CaseBResolution resolution,
                        std::size_t n_theta = st::kTrackAngleGridPoints,
                        std::size_t n_obs = st::kTrackObservedGridPoints) {
    const st::Grid En = st::energyGrid({});
    return {En, st::recoilAngleGrid(n_theta), observedGrid(n_obs), st::flatEnergyLogPrior(En),
            resolution};
}

// Standard normal CDF through erfc, accurate in both tails.
double normalCdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

struct Moments {
    double mean = 0.0;
    double std = 0.0;
};

std::vector<double> normalised(const std::vector<double>& log_post) {
    const double max_value = *std::max_element(log_post.begin(), log_post.end());
    std::vector<double> w(log_post.size());
    double sum = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) {
        w[i] = std::exp(log_post[i] - max_value);
        sum += w[i];
    }
    for (double& value : w) {
        value /= sum;
    }
    return w;
}

Moments gridMoments(const st::Grid& grid, const std::vector<double>& log_post) {
    const std::vector<double> w = normalised(log_post);
    double mean = 0.0;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        mean += w[i] * grid[i];
    }
    double var = 0.0;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        var += w[i] * (grid[i] - mean) * (grid[i] - mean);
    }
    return {mean, std::sqrt(var)};
}

// An acceptance histogram on 1 MeV x 5 deg bins over [10, 20] x [0, pi/2]
// whose ratio is shape(theta centre) in every E_n bin.
template <typename Shape> st::AcceptanceRatio acceptanceWithShape(Shape shape) {
    constexpr std::size_t kEnBins = 10;
    constexpr std::size_t kThetaBins = 18;
    st::AcceptanceHistogram h;
    for (std::size_t k = 0; k <= kEnBins; ++k) {
        h.En_edges.push_back(kEnMin + static_cast<double>(k));
    }
    for (std::size_t t = 0; t <= kThetaBins; ++t) {
        h.theta_edges.push_back(0.5 * pi * static_cast<double>(t) /
                                static_cast<double>(kThetaBins));
    }
    h.theta_edges.back() = 0.5 * pi;
    for (std::size_t k = 0; k < kEnBins; ++k) {
        const double En = 0.5 * (h.En_edges[k] + h.En_edges[k + 1]);
        for (std::size_t t = 0; t < kThetaBins; ++t) {
            const double lo = h.theta_edges[t];
            const double hi = h.theta_edges[t + 1];
            h.counts.push_back(shape(0.5 * (lo + hi)) * st::isotropicRecoilAngleMass(En, lo, hi));
        }
    }
    return st::AcceptanceRatio(h);
}

st::CaseBModel acceptanceModel(const st::AcceptanceRatio& ratio) {
    const st::Grid En = st::energyGrid({});
    const st::Grid theta = st::recoilAngleGrid(st::kTrackAngleGridPoints);
    return {En,
            theta,
            observedGrid(),
            st::flatEnergyLogPrior(En),
            kToyResolution,
            st::caseBAcceptanceLogWeights(ratio, En, theta)};
}

} // namespace

TEST_CASE("Case B: the track density of the kernel table is normalised, with the vMF first "
          "moment",
          "[stats][case_b]") {
    // Non-relativistic track density cos(theta_p)/pi convolved with the vMF:
    // E[cos theta_obs] = E[cos theta_p] (coth k - 1/k), E[cos theta_p] = 2/3
    // (addition theorem), as in the toy.
    const double sigma = 0.08;
    const double kappa = 1.0 / (sigma * sigma);
    const st::Grid obs = observedGrid(4001);
    const st::Grid theta = st::recoilAngleGrid(st::kTrackAngleGridPoints);
    const std::vector<double> table = st::logTrackKernelTable(obs, theta, sigma);
    REQUIRE(table.size() == obs.size() * theta.size());
    double norm = 0.0;
    double mean_cos = 0.0;
    for (std::size_t o = 0; o < obs.size(); ++o) {
        const double density =
            std::exp(st::logSumExp(table.data() + (o * theta.size()), theta.size()));
        const double weight = (o == 0 || o + 1 == obs.size()) ? 0.5 : 1.0;
        const double mass = weight * obs.step() * 2.0 * pi * std::sin(obs[o]) * density;
        norm += mass;
        mean_cos += mass * std::cos(obs[o]);
    }
    CHECK_THAT(norm, WithinAbs(1.0, 1e-4));
    CHECK_THAT(mean_cos, WithinAbs(2.0 / 3.0 * ((1.0 / std::tanh(kappa)) - (1.0 / kappa)), 1e-4));
}

TEST_CASE("Case B: sigma_theta -> 0 gives the sharp sum over E_n with the track density",
          "[stats][case_b]") {
    // Toy test: at sigma_theta = 0.003 the table at theta_obs equals the sum
    // over E_n at theta_p = theta_obs, with the relativistic track density,
    // where Ep / cos^2 stays inside the energy domain.
    const st::CaseBModel model = toyModel({0.10, 0.003}, 3000, 4001);
    const st::Grid& En = model.energyGrid();
    const std::vector<double> log_prior = st::flatEnergyLogPrior(En);
    std::size_t checked = 0;
    for (const double Ep_hat : {6.0, 8.0, 10.0}) {
        const std::vector<double> row = model.trackEnergyRow(Ep_hat);
        for (std::size_t t = 0; t < 40; ++t) {
            const double En_target = (kEnMin + (0.1 * (kEnMax - kEnMin))) +
                                     (static_cast<double>(t) / 39.0) * (0.6 * (kEnMax - kEnMin));
            const double theta_obs = std::acos(std::sqrt(std::min(Ep_hat / En_target, 1.0)));
            const double cos_obs = std::cos(theta_obs);
            const double En_implied = Ep_hat / (cos_obs * cos_obs);
            if (!(En_implied > kEnMin + (0.05 * (kEnMax - kEnMin)) && En_implied < 0.8 * kEnMax)) {
                continue;
            }
            std::vector<double> terms(En.size());
            for (std::size_t i = 0; i < En.size(); ++i) {
                const double z = (Ep_hat - kin::protonEnergy(En[i], theta_obs)) / 0.10;
                terms[i] = (-0.5 * z * z) + log_prior[i] + std::log(cos_obs / pi) +
                           st::logRecoilDensityRatio(theta_obs, En[i]);
            }
            const double sharp = st::logSumExp(terms.data(), terms.size());
            CHECK_THAT(model.interpolate(row.data(), theta_obs), WithinAbs(sharp, 2e-3));
            ++checked;
        }
    }
    CHECK(checked > 60);
}

TEST_CASE("Case B: Omega_n || z reduces to Case A up to the vMF kernel", "[stats][case_b]") {
    // Same theta_p prior on both sides: Case A rows carry the sin(2 theta_p)
    // cell mass times the relativistic factor, as the track density of Case B.
    // What is left is the 2D kernel against the 1D Gaussian, O(sigma_theta^2):
    // tolerances fixed before the run.
    const double sigma_theta = 0.02;
    const st::CaseBModel model_b = toyModel({0.10, sigma_theta});
    const st::Grid& En = model_b.energyGrid();
    const st::Grid& theta = model_b.recoilAngleGrid();
    const std::vector<double> mass = st::physicalRecoilAngleLogPrior(theta);
    std::vector<std::vector<double>> rows(En.size(), std::vector<double>(theta.size()));
    for (std::size_t i = 0; i < En.size(); ++i) {
        for (std::size_t j = 0; j < theta.size(); ++j) {
            rows[i][j] = mass[j] + st::logRecoilDensityRatio(theta[j], En[i]);
        }
    }
    const st::CaseAModel model_a(En, theta, rows, {0.10, sigma_theta});
    const double En_true = 15.0;
    for (const double theta_hat : {0.4, 0.7, 1.0}) {
        const double Ep_hat = kin::protonEnergy(En_true, theta_hat);
        const std::vector<double> post_a = model_a.logLikelihood(Ep_hat, theta_hat);
        const std::vector<double> post_b = model_b.energyLogPosterior(Ep_hat, theta_hat);
        const Moments a = gridMoments(En, post_a);
        const Moments b = gridMoments(En, post_b);
        CHECK(std::abs(a.mean - b.mean) <= 5.0 * sigma_theta * sigma_theta * a.mean);
        CHECK(std::abs(a.std - b.std) <= 0.05 * a.std);
        const std::vector<double> w_a = normalised(post_a);
        const std::vector<double> w_b = normalised(post_b);
        double total_variation = 0.0;
        for (std::size_t i = 0; i < En.size(); ++i) {
            total_variation += 0.5 * std::abs(w_a[i] - w_b[i]);
        }
        CHECK(total_variation <= 0.02);
    }
}

TEST_CASE("Case B: one event with a flat E_n prior has the analytic form of the toy",
          "[stats][case_b]") {
    // For sigma_theta -> 0, L ∝ cos(theta_p)/pi [Phi((20c - Ep)/s) -
    // Phi((10c - Ep)/s)] / c with c = cos^2(theta_p), the flat prior on E_n
    // integrated analytically. Within 30 degrees the smoothing by sigma_theta
    // and the relativistic factor stay below the tolerance of the toy; behind
    // pi/2 + 5 sigma_theta the likelihood is strongly suppressed, not -inf.
    const st::CaseBModel model = toyModel(kToyResolution);
    const std::vector<double> Ep_hat{11.0};
    const std::vector<geo::Direction> tracks{{0.0, 0.0, 1.0}};
    const std::vector<geo::Direction> candidates = st::sphereGrid(4000);
    const st::TrackEnergyTable table = model.trackEnergyTable(Ep_hat);
    const std::vector<double> loglik = model.summedLogLikelihood(table, tracks, candidates);
    REQUIRE(loglik.size() == candidates.size());
    double front_lo = std::numeric_limits<double>::infinity();
    double front_hi = -std::numeric_limits<double>::infinity();
    double front_max = -std::numeric_limits<double>::infinity();
    double rear_max = -std::numeric_limits<double>::infinity();
    std::size_t n_front = 0;
    std::size_t n_rear = 0;
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        REQUIRE(std::isfinite(loglik[c]));
        const double theta_p = geo::angleBetween(tracks[0], candidates[c]);
        if (theta_p < 30.0 * pi / 180.0) {
            const double c2 = std::cos(theta_p) * std::cos(theta_p);
            const double s = kToyResolution.sigma_Ep;
            const double analytic =
                std::log(std::cos(theta_p) / pi) - std::log(c2) +
                std::log(normalCdf((kEnMax * c2 - 11.0) / s) - normalCdf((kEnMin * c2 - 11.0) / s));
            front_lo = std::min(front_lo, loglik[c] - analytic);
            front_hi = std::max(front_hi, loglik[c] - analytic);
            front_max = std::max(front_max, loglik[c]);
            ++n_front;
        } else if (theta_p > (pi / 2.0) + (5.0 * kToyResolution.sigma_theta)) {
            rear_max = std::max(rear_max, loglik[c]);
            ++n_rear;
        }
    }
    REQUIRE(n_front > 100);
    REQUIRE(n_rear > 100);
    CHECK(front_hi - front_lo < 0.02);
    CHECK(rear_max < front_max - 20.0);
}

TEST_CASE("Case B: summed log-likelihood is the sum of the interpolated rows", "[stats][case_b]") {
    const st::CaseBModel model = toyModel(kToyResolution);
    const std::vector<double> Ep_hat{4.0, 11.0, -0.05};
    const std::vector<geo::Direction> tracks{
        {0.0, 0.0, 1.0}, geo::normalized({1.0, 2.0, 0.5}), geo::normalized({-0.3, 0.1, 0.9})};
    const std::vector<geo::Direction> candidates = st::sphereGrid(50);
    const st::TrackEnergyTable table = model.trackEnergyTable(Ep_hat);
    const std::vector<double> summed = model.summedLogLikelihood(table, tracks, candidates);
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        double expected = 0.0;
        for (std::size_t k = 0; k < tracks.size(); ++k) {
            const std::vector<double> row = model.trackEnergyRow(Ep_hat[k]);
            expected += model.interpolate(row.data(), geo::angleBetween(tracks[k], candidates[c]));
        }
        CHECK_THAT(summed[c], WithinAbs(expected, 0.0));
    }
    // At a grid node the interpolation returns the node value.
    const std::vector<double> row = model.trackEnergyRow(11.0);
    const st::Grid& obs = model.observedAngleGrid();
    CHECK_THAT(model.interpolate(row.data(), obs[17]), WithinAbs(row[17], 1e-12));
    CHECK_THAT(model.interpolate(row.data(), pi), WithinAbs(row.back(), 1e-6));
}

TEST_CASE("Case B: the E_n window leaves the table unchanged", "[stats][case_b]") {
    // The window skips only terms below L - C; the full sum differs by at
    // most n_En exp(-50) relative, far below 1e-12 in log.
    const st::CaseBModel model = toyModel(kToyResolution);
    const double Ep_hat = 9.0;
    const std::vector<double> row = model.trackEnergyRow(Ep_hat);
    const st::Grid& obs = model.observedAngleGrid();
    for (const std::size_t o :
         {std::size_t{0}, std::size_t{200}, std::size_t{350}, std::size_t{600}}) {
        const std::vector<double> post = model.energyLogPosterior(Ep_hat, obs[o]);
        // Both sides carry the E_n prior: the row through the weights, the
        // diagnostic explicitly.
        const double direct = st::logSumExp(post.data(), post.size());
        CHECK_THAT(row[o], WithinAbs(direct, 1e-9));
    }
}

TEST_CASE("Case B: invalid inputs throw", "[stats][case_b]") {
    const st::Grid En = st::energyGrid({});
    const st::Grid theta = st::recoilAngleGrid(st::kTrackAngleGridPoints);
    const std::vector<double> prior = st::flatEnergyLogPrior(En);
    CHECK_THROWS_AS(st::CaseBModel(En, theta, observedGrid(), prior, {0.0, 0.08}),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::CaseBModel(En, theta, observedGrid(), prior, {0.1, -1.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(
        st::CaseBModel(En, st::Grid(0.0, 2.0, 100), observedGrid(), prior, kToyResolution),
        std::invalid_argument);
    CHECK_THROWS_AS(st::CaseBModel(En, theta, st::Grid(0.0, 3.0, 801), prior, kToyResolution),
                    std::invalid_argument);
    CHECK_THROWS_AS(
        st::CaseBModel(st::Grid(0.0, 20.0, 500), theta, observedGrid(), prior, kToyResolution),
        std::invalid_argument);
    CHECK_THROWS_AS(
        st::CaseBModel(En, theta, observedGrid(), std::vector<double>(3, 0.0), kToyResolution),
        std::invalid_argument);
    std::vector<double> bad = prior;
    bad[4] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(st::CaseBModel(En, theta, observedGrid(), bad, kToyResolution),
                    std::invalid_argument);
    const std::vector<double> empty(En.size(), -std::numeric_limits<double>::infinity());
    CHECK_THROWS_AS(st::CaseBModel(En, theta, observedGrid(), empty, kToyResolution),
                    std::invalid_argument);
    const st::CaseBModel model(En, theta, observedGrid(), prior, kToyResolution);
    CHECK_THROWS_AS(model.trackEnergyRow(std::numeric_limits<double>::infinity()),
                    std::invalid_argument);
    const st::TrackEnergyTable table = model.trackEnergyTable({11.0});
    CHECK_THROWS_AS(model.summedLogLikelihood(table, {}, st::sphereGrid(10)),
                    std::invalid_argument);
}

TEST_CASE("Case B: an acceptance r = 1 reproduces the ideal model", "[stats][case_b]") {
    const st::CaseBModel ideal = toyModel(kToyResolution);
    const st::CaseBModel flat = acceptanceModel(acceptanceWithShape([](double) { return 1.0; }));
    CHECK_FALSE(ideal.hasAcceptance());
    CHECK(flat.hasAcceptance());
    const std::vector<double> Ep_hat{4.0, 11.0, 17.5, -0.05};
    const std::vector<geo::Direction> tracks{{0.0, 0.0, 1.0},
                                             geo::normalized({1.0, 2.0, 0.5}),
                                             geo::normalized({-0.3, 0.1, 0.9}),
                                             geo::normalized({0.2, -0.7, 0.1})};
    const std::vector<geo::Direction> candidates = st::sphereGrid(200);
    const std::vector<double> a =
        ideal.summedLogLikelihood(ideal.trackEnergyTable(Ep_hat), tracks, candidates);
    const std::vector<double> b =
        flat.summedLogLikelihood(flat.trackEnergyTable(Ep_hat), tracks, candidates);
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        CHECK_THAT(b[c], WithinAbs(a[c], 1e-9));
    }
    const std::vector<double> pa = ideal.energyLogPosterior(11.0, 0.4);
    const std::vector<double> pb = flat.energyLogPosterior(11.0, 0.4);
    for (std::size_t i = 0; i < pa.size(); ++i) {
        CHECK_THAT(pb[i], WithinAbs(pa[i], 1e-9));
    }
}

TEST_CASE("Case B: a forward acceptance deficit moves the E_n posterior up", "[stats][case_b]") {
    // At fixed Ep_hat a small recoil angle goes with a small E_n (Ep = E_n
    // cos^2 theta_p in the non-relativistic limit): fewer forward recoils
    // shift the posterior of E_n to larger values.
    constexpr double kForwardEdge = 0.5;
    const st::CaseBModel ideal = toyModel(kToyResolution);
    const st::CaseBModel step = acceptanceModel(
        acceptanceWithShape([](double theta) { return theta < kForwardEdge ? 0.3 : 1.0; }));
    const st::Grid& En = ideal.energyGrid();
    const double Ep_hat = 11.0;
    const double theta_obs = 0.5;
    const Moments a = gridMoments(En, ideal.energyLogPosterior(Ep_hat, theta_obs));
    const Moments b = gridMoments(En, step.energyLogPosterior(Ep_hat, theta_obs));
    CHECK(b.mean > a.mean + 0.1);
    // The table rows follow the same diagnostic, as for the ideal model.
    const std::vector<double> row = step.trackEnergyRow(Ep_hat);
    const st::Grid& obs = step.observedAngleGrid();
    const std::vector<double> post = step.energyLogPosterior(Ep_hat, obs[300]);
    CHECK_THAT(row[300], WithinAbs(st::logSumExp(post.data(), post.size()), 1e-9));
}

TEST_CASE("Case B: an invalid acceptance table throws", "[stats][case_b]") {
    const st::Grid En = st::energyGrid({});
    const st::Grid theta = st::recoilAngleGrid(st::kTrackAngleGridPoints);
    const std::vector<double> prior = st::flatEnergyLogPrior(En);
    CHECK_THROWS_AS(st::CaseBModel(En, theta, observedGrid(), prior, kToyResolution,
                                   std::vector<double>(7, 0.0)),
                    std::invalid_argument);
    std::vector<double> bad(theta.size() * En.size(), 0.0);
    bad[12] = -std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(st::CaseBModel(En, theta, observedGrid(), prior, kToyResolution, bad),
                    std::invalid_argument);
}
