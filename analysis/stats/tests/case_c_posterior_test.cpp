#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_b_posterior.hpp"
#include "riptide/stats/case_c.hpp"
#include "riptide/stats/case_c_posterior.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/log_sum_exp.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"

namespace st = riptide::stats;
namespace geo = riptide::geometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr double kMinusInf = -std::numeric_limits<double>::infinity();

st::HyperGrid smallGrid() {
    return {st::Grid{12.0, 16.0, 5}, st::Grid{std::log(0.5), std::log(2.0), 4}};
}

// All the mass on one point of the grid.
st::HyperPosterior pointMass(std::size_t index) {
    st::HyperPosterior posterior{smallGrid(), {}};
    posterior.log_posterior.assign(posterior.grid.size(), kMinusInf);
    posterior.log_posterior[index] = 0.0;
    return posterior;
}

double massOf(const std::vector<double>& log_mass) {
    return std::exp(st::logSumExp(log_mass));
}

} // namespace

TEST_CASE("predictive moments of a point mass are its mu_E and sigma_E", "[case_c]") {
    const st::HyperPosterior posterior = pointMass(13);
    const st::PredictiveEnergy predictive = st::predictiveEnergyMoments(posterior);
    CHECK_THAT(predictive.mean, WithinRel(posterior.grid.muAt(13), 1e-15));
    CHECK_THAT(predictive.std, WithinRel(posterior.grid.sigmaAt(13), 1e-12));
}

TEST_CASE("predictive variance adds the spread of mu_E to sigma_E^2", "[case_c]") {
    // Equal mass on (mu[0], sigma[1]) and (mu[4], sigma[1]): Var = sigma^2 +
    // ((mu[4] - mu[0]) / 2)^2.
    st::HyperPosterior posterior = pointMass(1);
    const std::size_t other = posterior.grid.index(4, 1);
    posterior.log_posterior[other] = 0.0;
    const double half_gap = 0.5 * (posterior.grid.mu[4] - posterior.grid.mu[0]);
    const double sigma = posterior.grid.sigmaAt(1);
    const st::PredictiveEnergy predictive = st::predictiveEnergyMoments(posterior);
    CHECK_THAT(predictive.mean, WithinRel(14.0, 1e-15));
    CHECK_THAT(predictive.std,
               WithinRel(std::sqrt((sigma * sigma) + (half_gap * half_gap)), 1e-12));
}

TEST_CASE("the marginals of a separable posterior are its factors", "[case_c]") {
    st::HyperPosterior posterior{smallGrid(), {}};
    const st::HyperGrid& grid = posterior.grid;
    const std::array<double, 5> f{-1.0, 0.3, 2.0, -0.4, -3.0};
    const std::array<double, 4> g{0.5, -2.0, 1.0, 0.0};
    posterior.log_posterior.resize(grid.size());
    for (std::size_t i = 0; i < f.size(); ++i) {
        for (std::size_t j = 0; j < g.size(); ++j) {
            posterior.log_posterior[grid.index(i, j)] = f[i] + g[j] + 7.0;
        }
    }
    const std::vector<double> mu = st::muLogMarginal(posterior);
    const std::vector<double> ls = st::logSigmaLogMarginal(posterior);
    REQUIRE(mu.size() == f.size());
    REQUIRE(ls.size() == g.size());
    CHECK_THAT(massOf(mu), WithinRel(1.0, 1e-14));
    CHECK_THAT(massOf(ls), WithinRel(1.0, 1e-14));
    const double log_f = st::logSumExp(std::vector<double>(f.begin(), f.end()));
    const double log_g = st::logSumExp(std::vector<double>(g.begin(), g.end()));
    for (std::size_t i = 0; i < f.size(); ++i) {
        CHECK_THAT(mu[i], WithinAbs(f[i] - log_f, 1e-13));
    }
    for (std::size_t j = 0; j < g.size(); ++j) {
        CHECK_THAT(ls[j], WithinAbs(g[j] - log_g, 1e-13));
    }
}

TEST_CASE("a symmetric marginal has its median at the centre, inside every interval", "[case_c]") {
    const st::Grid axis{-2.0, 2.0, 41};
    std::vector<double> log_marginal(axis.size());
    for (std::size_t i = 0; i < axis.size(); ++i) {
        log_marginal[i] = -0.5 * axis[i] * axis[i] / 0.25;
    }
    const st::HyperparameterSummary summary = st::summarizeHyperparameter(axis, log_marginal);
    CHECK_THAT(summary.median, WithinAbs(0.0, 1e-12));
    CHECK_THAT(summary.mean, WithinAbs(0.0, 1e-12));
    for (const st::Interval& interval : summary.intervals) {
        CHECK(interval.contains(summary.median));
        CHECK_THAT(interval.lo, WithinAbs(-interval.hi, 1e-12));
    }
    CHECK(summary.intervals[0].hi < summary.intervals[1].hi);
    CHECK(summary.intervals[1].hi < summary.intervals[2].hi);
}

TEST_CASE("summaries reject a log-posterior off its grid", "[case_c]") {
    st::HyperPosterior posterior = pointMass(0);
    posterior.log_posterior.pop_back();
    CHECK_THROWS_AS(st::muLogMarginal(posterior), std::invalid_argument);
    CHECK_THROWS_AS(st::predictiveEnergyMoments(posterior), std::invalid_argument);
}

TEST_CASE("the two-stage analysis reproduces the toy model on a fixed experiment", "[case_c]") {
    // 50 events of the Case B generator with E_n uniform on [13, 16] MeV.
    // Reference values: the toy-model refine_hyperparameters (both hyperpriors),
    // posterior_mean_std, credible_interval (level 0 for the median) and
    // predictive_energy_moments on these events and the omega_0 of the C++
    // stage 1; refine_shared_direction and its hierarchical form give the
    // same pixel as here for omega_0 and omega_hat.
    constexpr double kPi = 3.14159265358979323846;
    constexpr std::size_t kEvents = 50;
    constexpr std::uint64_t kSeed = 20260907;
    const st::CaseBResolution resolution{0.10, 0.08};
    const st::CaseCModel model(st::energyGrid({}), st::recoilAngleGrid(st::kTrackAngleGridPoints),
                               st::Grid{0.0, kPi, st::kTrackObservedGridPoints}, resolution);
    const geo::Direction omega_n = geo::normalized({0.3, -0.5, 0.8});
    const st::CaseBPseudoSample sample = st::caseBPseudoSample(
        st::EnergyDomain{13.0, 16.0, st::kEnergyGridPoints}, omega_n, kEvents, resolution, kSeed);
    const st::HyperLogPrior baseline = st::uniformLogSigmaHyperPrior;
    const st::HyperLogPrior alternate = st::uniformSigmaHyperPrior;
    const st::CaseCAnalysis analysis = st::analyseCaseCExperiment(
        model, sample.reco_Ep, sample.reco_track, st::sphereGrid(st::kSphereGridPixels),
        st::uniformDirectionLogPrior, baseline, &alternate);

    constexpr double kTolerance = 1e-3;
    const auto check = [&](const st::HyperparameterSummary& got, const std::array<double, 9>& toy) {
        CHECK_THAT(got.mean, WithinRel(toy[0], kTolerance));
        CHECK_THAT(got.std, WithinRel(toy[1], kTolerance));
        CHECK_THAT(got.median, WithinRel(toy[2], kTolerance));
        for (std::size_t l = 0; l < got.intervals.size(); ++l) {
            CHECK_THAT(got.intervals[l].lo, WithinRel(toy[3 + (2 * l)], kTolerance));
            CHECK_THAT(got.intervals[l].hi, WithinRel(toy[4 + (2 * l)], kTolerance));
        }
    };
    check(analysis.hyper.mu, {14.11973252, 0.2660021854, 14.1020659, 13.85933918, 14.38372353,
                              13.6890519, 14.60156525, 13.61223884, 14.73068297});
    check(analysis.hyper.log_sigma,
          {-1.014835064, 1.354696074, -0.5011811429, -2.625396835, 0.08977398046, -4.043748634,
           0.384033792, -4.369869877, 0.5119930886});
    REQUIRE(analysis.alternate.has_value());
    check(analysis.alternate->mu, {14.16315773, 0.5109795036, 14.15019255, 13.82357974, 14.5132476,
                                   13.6011589, 14.79086349, 13.42148964, 14.94909953});
    check(analysis.alternate->log_sigma,
          {-0.09373456643, 0.7378804963, -0.06505474596, -0.5933555213, 0.3695434171, -1.157393953,
           0.7118711444, -1.636301693, 1.055228489});
    CHECK_THAT(analysis.predictive.mean, WithinRel(14.119732520728105, kTolerance));
    CHECK_THAT(analysis.predictive.std, WithinRel(0.8416323338277383, kTolerance));

    // Both maps are pixels of the lattices, so a tolerance far below one
    // cap pixel (~1e-3 rad here) checks that the pixel is the same.
    constexpr double kSamePixel = 1e-6;
    CHECK(geo::angleBetween(analysis.omega_0, {0.33534217335322641, -0.48538072203872773,
                                               0.8074349394495498}) < kSamePixel);
    CHECK(geo::angleBetween(analysis.direction.omega_hat,
                            {0.32929233466321128, -0.46721085113050437, 0.82053676268523246}) <
          kSamePixel);
}

TEST_CASE("the analysis without an alternate prior has no alternate summary", "[case_c]") {
    constexpr double kPi = 3.14159265358979323846;
    const st::CaseBResolution resolution{0.10, 0.08};
    const st::CaseCModel model(st::energyGrid({}), st::recoilAngleGrid(st::kTrackAngleGridPoints),
                               st::Grid{0.0, kPi, st::kTrackObservedGridPoints}, resolution);
    const geo::Direction omega_n{0.0, 0.0, 1.0};
    const st::CaseBPseudoSample sample = st::caseBPseudoSample(
        st::EnergyDomain{13.0, 16.0, st::kEnergyGridPoints}, omega_n, 20, resolution, 20260907);
    const st::CaseCAnalysis analysis = st::analyseCaseCExperiment(
        model, sample.reco_Ep, sample.reco_track, st::sphereGrid(st::kSphereGridPixels),
        st::uniformDirectionLogPrior, st::uniformLogSigmaHyperPrior);
    CHECK_FALSE(analysis.alternate.has_value());
    CHECK(analysis.predictive.std > 0.0);
    CHECK(geo::angleBetween(analysis.direction.omega_hat, omega_n) < 0.5);
}
