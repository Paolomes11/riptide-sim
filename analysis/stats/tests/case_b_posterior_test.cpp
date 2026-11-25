#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_b_posterior.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/log_sum_exp.hpp"
#include "riptide/stats/priors.hpp"

namespace st = riptide::stats;
namespace kin = riptide::kinematics;
namespace geo = riptide::geometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;
using riptide::units::pi;

namespace {

constexpr st::CaseBResolution kToyResolution{0.10, 0.08};

st::CaseBModel toyModel() {
    const st::Grid En = st::energyGrid({});
    return {En, st::recoilAngleGrid(st::kTrackAngleGridPoints),
            st::Grid{0.0, pi, st::kTrackObservedGridPoints}, st::flatEnergyLogPrior(En),
            kToyResolution};
}

const geo::Direction kTruth = geo::normalized({1.0, 2.0, 2.0});

// Noise-free sample symmetric about kTruth: three recoil angles, eight
// azimuths each, Ep exactly the kinematic value at 15 MeV. By symmetry the
// likelihood peaks on kTruth itself, so the MAP is off only by the grid.
struct Sample {
    std::vector<double> Ep_hat;
    std::vector<geo::Direction> tracks;
};

Sample symmetricSample() {
    constexpr double kEn = 15.0;
    constexpr std::size_t kAzimuths = 8;
    Sample sample;
    for (const double theta_p : {0.3, 0.6, 0.9}) {
        for (std::size_t a = 0; a < kAzimuths; ++a) {
            const double phi = 2.0 * pi * static_cast<double>(a) / static_cast<double>(kAzimuths);
            sample.tracks.push_back(geo::directionAroundAxis(kTruth, theta_p, phi));
            sample.Ep_hat.push_back(kin::protonEnergy(kEn, theta_p));
        }
    }
    return sample;
}

// Hand-built posterior on a cap about +z: a 2D Gaussian of width s per
// tangent component, log p = -alpha^2 / (2 s^2), alpha the angle from +z.
constexpr double kCapRadius = 0.3;
constexpr double kWidth = 0.03;
const geo::Direction kZ{0.0, 0.0, 1.0};

st::DirectionPosterior gaussianCapPosterior() {
    st::DirectionPosterior posterior;
    posterior.cap_radius = kCapRadius;
    posterior.cap = st::directionCapGrid(kZ, kCapRadius, st::kCapGridPixels);
    posterior.cap_log_posterior.resize(posterior.cap.size());
    for (std::size_t i = 0; i < posterior.cap.size(); ++i) {
        const double alpha = geo::angleBetween(posterior.cap[i], kZ);
        posterior.cap_log_posterior[i] = -alpha * alpha / (2.0 * kWidth * kWidth);
    }
    // Pixel 0 is the nearest to the centre (cos alpha_0 is the largest).
    posterior.map_index = 0;
    posterior.omega_hat = posterior.cap[0];
    return posterior;
}

// The direction at angle `alpha` from +z at which the Gaussian HPD region
// has mass m strictly inside: 1 - exp(-alpha^2 / (2 s^2)) = m.
double radiusOfMass(double m) {
    return kWidth * std::sqrt(-2.0 * std::log(1.0 - m));
}

} // namespace

TEST_CASE("uniform direction prior is -log n and sums to one", "[stats][case_b_posterior]") {
    const auto grid = st::sphereGrid(st::kSphereGridPixels);
    const auto log_w = st::uniformDirectionLogPrior(grid);
    REQUIRE(log_w.size() == grid.size());
    for (const double value : log_w) {
        REQUIRE_THAT(value, WithinAbs(-std::log(3000.0), 1e-15));
    }
    REQUIRE_THAT(st::logSumExp(log_w), WithinAbs(0.0, 1e-12));
    REQUIRE_THROWS_AS(st::uniformDirectionLogPrior({}), std::invalid_argument);
}

TEST_CASE("vMF direction prior: normalised, uniform at kappa 0, mean cosine coth k - 1/k",
          "[stats][case_b_posterior]") {
    // Declared before the run: normalisation and the kappa = 0 limit to 1e-12;
    // the first moment on a 20000-pixel equal-area lattice to 1e-3 (the
    // quadrature error of a smooth function on that lattice is far smaller).
    const auto grid = st::sphereGrid(20000);
    const geo::Direction axis = geo::normalized({0.3, -0.5, 0.8});

    const auto uniform = st::uniformDirectionLogPrior(grid);
    const auto at_zero = st::vmfDirectionLogPrior(grid, axis, 0.0);
    for (std::size_t i = 0; i < grid.size(); ++i) {
        REQUIRE_THAT(at_zero[i], WithinAbs(uniform[i], 1e-12));
    }

    for (const double kappa : {1.0, 5.0}) {
        const auto log_w = st::vmfDirectionLogPrior(grid, axis, kappa);
        REQUIRE_THAT(st::logSumExp(log_w), WithinAbs(0.0, 1e-12));
        double mean_cos = 0.0;
        for (std::size_t i = 0; i < grid.size(); ++i) {
            mean_cos += std::exp(log_w[i]) * geo::dot(grid[i], axis);
        }
        const double expected = 1.0 / std::tanh(kappa) - 1.0 / kappa;
        REQUIRE_THAT(mean_cos, WithinAbs(expected, 1e-3));
    }

    REQUIRE_THROWS_AS(st::vmfDirectionLogPrior(grid, axis, -1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(
        st::vmfDirectionLogPrior(grid, axis, std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);
    REQUIRE_THROWS_AS(st::vmfDirectionLogPrior({}, axis, 1.0), std::invalid_argument);
}

TEST_CASE("direction log-posterior check", "[stats][case_b_posterior]") {
    constexpr double kInf = std::numeric_limits<double>::infinity();
    REQUIRE_NOTHROW(st::checkDirectionLogPosterior({-kInf, -3.0, -kInf}));
    REQUIRE_THROWS_AS(st::checkDirectionLogPosterior({-kInf, -kInf}), std::domain_error);
    REQUIRE_THROWS_AS(st::checkDirectionLogPosterior({}), std::domain_error);
    REQUIRE_THROWS_AS(
        st::checkDirectionLogPosterior({0.0, std::numeric_limits<double>::quiet_NaN()}),
        std::domain_error);
    REQUIRE_THROWS_AS(st::checkDirectionLogPosterior({0.0, kInf}), std::domain_error);
}

TEST_CASE("summaries on a Gaussian posterior with known width", "[stats][case_b_posterior]") {
    // Declared before the run. sigma: sqrt(E[alpha^2]) = sqrt(2) s for a 2D
    // Gaussian, within 2% (the cap edge at 10 s truncates nothing visible;
    // the MAP pixel sits ~2e-3 rad off the centre). Coverage: truths at HPD
    // mass 0.50 / 0.80 / 0.93 / 0.97, at least 0.02 from every level, so
    // the discretisation of the mass (~1e-3) cannot flip an outcome.
    const st::DirectionPosterior posterior = gaussianCapPosterior();

    struct Case {
        double mass;
        std::array<bool, 3> covered;
    };
    for (const Case& c : {Case{0.50, {true, true, true}}, Case{0.80, {false, true, true}},
                          Case{0.93, {false, false, true}}, Case{0.97, {false, false, false}}}) {
        CAPTURE(c.mass);
        const geo::Direction truth = geo::directionAroundAxis(kZ, radiusOfMass(c.mass), 1.0);
        const st::DirectionSummary s = st::summarizeDirection(posterior, truth);
        REQUIRE_THAT(s.sigma, WithinRel(std::sqrt(2.0) * kWidth, 0.02));
        REQUIRE_THAT(s.error, WithinAbs(geo::angleBetween(truth, posterior.omega_hat), 0.0));
        REQUIRE_THAT(s.pull, WithinRel(s.error / s.sigma, 1e-15));
        REQUIRE(s.truth_in_cap);
        for (std::size_t l = 0; l < 3; ++l) {
            CAPTURE(l);
            REQUIRE(s.covered[l] == c.covered[l]);
        }
    }

    // A truth well outside the cap is never covered.
    const geo::Direction outside = geo::directionAroundAxis(kZ, 0.5, 0.0);
    const st::DirectionSummary s = st::summarizeDirection(posterior, outside);
    REQUIRE_FALSE(s.truth_in_cap);
    for (const bool covered : s.covered) {
        REQUIRE_FALSE(covered);
    }
}

TEST_CASE("Case B posterior: finite on the sphere, MAP at the truth, reproducible",
          "[stats][case_b_posterior]") {
    // Declared before the run. The likelihood is finite on every coarse
    // pixel (the vMF kernel has no zero). With the symmetric noise-free
    // sample the MAP is within 0.015 rad of the truth (cap pixels are
    // ~6e-3 rad for a radius ~0.3). Repeating the call gives the same bits.
    const st::CaseBModel model = toyModel();
    const Sample sample = symmetricSample();
    const st::TrackEnergyTable table = model.trackEnergyTable(sample.Ep_hat);
    const auto coarse_grid = st::sphereGrid(st::kSphereGridPixels);

    const auto coarse = model.summedLogLikelihood(table, sample.tracks, coarse_grid);
    for (const double value : coarse) {
        REQUIRE(std::isfinite(value));
    }

    const st::DirectionPosterior posterior = st::refineSharedDirection(
        model, table, sample.tracks, coarse_grid, st::uniformDirectionLogPrior);
    REQUIRE(posterior.cap.size() == st::kCapGridPixels);
    REQUIRE(posterior.cap_radius > 0.0);
    REQUIRE(posterior.cap_radius <= pi);
    REQUIRE(posterior.map_index < posterior.cap.size());
    for (const double value : posterior.cap_log_posterior) {
        REQUIRE(value <= posterior.cap_log_posterior[posterior.map_index]);
    }
    REQUIRE(geo::angleBetween(posterior.omega_hat, kTruth) < 0.015);

    const st::DirectionSummary summary = st::summarizeDirection(posterior, kTruth);
    REQUIRE(summary.truth_in_cap);
    REQUIRE(summary.covered[0]);

    const st::DirectionPosterior again = st::refineSharedDirection(
        model, table, sample.tracks, coarse_grid, st::uniformDirectionLogPrior);
    REQUIRE(again.map_index == posterior.map_index);
    REQUIRE_THAT(again.cap_radius, WithinULP(posterior.cap_radius, 0));
    for (std::size_t i = 0; i < posterior.cap.size(); ++i) {
        REQUIRE_THAT(again.cap_log_posterior[i], WithinULP(posterior.cap_log_posterior[i], 0));
    }
}

TEST_CASE("Case B posterior: a sharp vMF prior pulls the MAP, bad inputs throw",
          "[stats][case_b_posterior]") {
    // Declared before the run: a prior of width 1/sqrt(2000) ~ 0.022 rad on an
    // axis 0.1 rad from the truth is narrower than the likelihood, so the MAP
    // moves toward the axis: closer to it than the truth is, and more than
    // 0.015 rad (the uniform-prior bound above) from the truth.
    const st::CaseBModel model = toyModel();
    const Sample sample = symmetricSample();
    const st::TrackEnergyTable table = model.trackEnergyTable(sample.Ep_hat);
    const auto coarse_grid = st::sphereGrid(st::kSphereGridPixels);

    const geo::Direction axis = geo::directionAroundAxis(kTruth, 0.1, 0.7);
    const st::DirectionLogPrior vmf = [&axis](const std::vector<geo::Direction>& grid) {
        return st::vmfDirectionLogPrior(grid, axis, 2000.0);
    };
    const st::DirectionPosterior posterior =
        st::refineSharedDirection(model, table, sample.tracks, coarse_grid, vmf);
    REQUIRE(geo::angleBetween(posterior.omega_hat, axis) < geo::angleBetween(kTruth, axis));
    REQUIRE(geo::angleBetween(posterior.omega_hat, kTruth) > 0.015);

    const st::DirectionLogPrior wrong_size = [](const std::vector<geo::Direction>& grid) {
        return std::vector<double>(grid.size() + 1, 0.0);
    };
    REQUIRE_THROWS_AS(
        st::refineSharedDirection(model, table, sample.tracks, coarse_grid, wrong_size),
        std::invalid_argument);

    std::vector<geo::Direction> fewer_tracks = sample.tracks;
    fewer_tracks.pop_back();
    REQUIRE_THROWS_AS(st::refineSharedDirection(model, table, fewer_tracks, coarse_grid,
                                                st::uniformDirectionLogPrior),
                      std::invalid_argument);
}
