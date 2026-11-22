#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/combine.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/log_sum_exp.hpp"
#include "riptide/stats/priors.hpp"

namespace st = riptide::stats;
namespace kin = riptide::kinematics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace {

// Resolutions of the toy reference cases (D31).
constexpr st::CaseAResolution kToyResolution{0.10, 0.08};

struct Moments {
    double mean = 0.0;
    double std = 0.0;
};

// Weighted moments on the grid, as the toy reference tests compute them.
Moments gridMoments(const st::Grid& grid, const std::vector<double>& log_post) {
    const double max_value = *std::max_element(log_post.begin(), log_post.end());
    double w_sum = 0.0;
    double w_x = 0.0;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        const double w = std::exp(log_post[i] - max_value);
        w_sum += w;
        w_x += w * grid[i];
    }
    const double mean = w_x / w_sum;
    double w_dx2 = 0.0;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        const double dx = grid[i] - mean;
        w_dx2 += std::exp(log_post[i] - max_value) * dx * dx;
    }
    return {mean, std::sqrt(w_dx2 / w_sum)};
}

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

std::size_t argmax(const std::vector<double>& values) {
    return static_cast<std::size_t>(
        std::distance(values.begin(), std::max_element(values.begin(), values.end())));
}

st::CaseAModel flatThetaModel(const st::Grid& En_grid, st::CaseAResolution resolution) {
    const st::Grid theta = st::recoilAngleGrid();
    return {En_grid, theta, st::flatRecoilAngleLogPrior(theta), resolution};
}

st::CaseAModel physicalModel(const st::Grid& En_grid, st::CaseAResolution resolution) {
    const st::Grid theta = st::recoilAngleGrid();
    return {En_grid, theta, st::physicalRecoilAngleLogPrior(theta), resolution};
}

} // namespace

TEST_CASE("R1: one event gives E_n = 1.82 +/- 0.21 MeV", "[stats][case_a][oracle]") {
    // The low-energy grid of the toy reference test: R1 checks the engine,
    // not the physical domain, which does not contain 1.82 MeV.
    const st::Grid En(0.5, 6.0, 500);
    const std::vector<double> log_prior_En = st::flatEnergyLogPrior(En);

    // R1 is defined by linearisation, without the physical theta_p prior.
    const st::CaseAModel flat = flatThetaModel(En, kToyResolution);
    const std::vector<double> post = st::singleEventLogPosterior(flat, 1.40, 0.50, log_prior_En);
    CHECK_THAT(En[argmax(post)], WithinAbs(1.82, 0.03));
    CHECK_THAT(gridMoments(En, post).std, WithinAbs(0.21, 0.01));

    // With sin(2 theta_p) the toy finds sigma = 0.222 MeV.
    const st::CaseAModel physical = physicalModel(En, kToyResolution);
    const std::vector<double> post_physical =
        st::singleEventLogPosterior(physical, 1.40, 0.50, log_prior_En);
    CHECK_THAT(gridMoments(En, post_physical).std, WithinAbs(0.222, 0.005));
}

TEST_CASE("R2: two events combine to E_n = 2.40 +/- 0.18 MeV", "[stats][combine][oracle]") {
    const st::Grid En(0.5, 6.0, 2000);
    st::LogLikelihoodMatrix ll;
    ll.events = 2;
    ll.points = En.size();
    ll.values.resize(ll.events * ll.points);
    for (std::size_t i = 0; i < En.size(); ++i) {
        const double z1 = (En[i] - 2.70) / 0.45;
        const double z2 = (En[i] - 2.35) / 0.20;
        ll.values[i] = -0.5 * z1 * z1;
        ll.values[ll.points + i] = -0.5 * z2 * z2;
    }
    const std::vector<double> combined =
        st::combineLogLikelihoods(ll, std::vector<double>(En.size(), 0.0));
    CHECK_THAT(En[argmax(combined)], WithinAbs(2.40, 0.03));
    CHECK_THAT(gridMoments(En, combined).std, WithinAbs(0.18, 0.01));
}

TEST_CASE("a noise-free event at 15 MeV is recovered within 3 sigma", "[stats][case_a]") {
    constexpr double kEn = 15.0;
    constexpr double kTheta = 0.4;
    const st::Grid En = st::energyGrid(st::EnergyDomain{});
    const st::CaseAModel model = flatThetaModel(En, kToyResolution);
    const std::vector<double> post = st::singleEventLogPosterior(
        model, kin::protonEnergy(kEn, kTheta), kTheta, st::flatEnergyLogPrior(En));

    // Linearised propagation of Ep = En cos^2(theta).
    const double cos2 = std::cos(kTheta) * std::cos(kTheta);
    const double sigma_En = std::hypot(kToyResolution.sigma_Ep / cos2,
                                       2.0 * kEn * std::tan(kTheta) * kToyResolution.sigma_theta);
    CHECK_THAT(En[argmax(post)], WithinAbs(kEn, 3.0 * sigma_En));
}

TEST_CASE("Case A limits", "[stats][case_a]") {
    const st::Grid En = st::energyGrid(st::EnergyDomain{});
    const std::vector<double> log_prior_En = st::flatEnergyLogPrior(En);

    SECTION("noise-free data with a tiny resolution pick the true cell") {
        const st::CaseAModel model = physicalModel(En, st::CaseAResolution{1e-4, 1e-5});
        const std::size_t i_true = 137;
        const std::size_t j_true = 211;
        const double theta = model.recoilAngleGrid()[j_true];
        const std::vector<double> post = st::singleEventLogPosterior(
            model, kin::protonEnergy(En[i_true], theta), theta, log_prior_En);
        CHECK(argmax(post) == i_true);
        CHECK(normalised(post)[i_true] > 0.99);
    }

    SECTION("a huge resolution gives back the prior") {
        const st::CaseAModel model = physicalModel(En, st::CaseAResolution{1e3, 1e3});
        const std::vector<double> post =
            normalised(st::singleEventLogPosterior(model, 7.0, 0.6, log_prior_En));
        const std::vector<double> prior = normalised(log_prior_En);
        for (std::size_t i = 0; i < En.size(); ++i) {
            CHECK_THAT(post[i], WithinRel(prior[i], 1e-3));
        }
    }

    SECTION("N identical events weigh N times one event") {
        constexpr std::size_t kEvents = 25;
        const st::CaseAModel model = physicalModel(En, kToyResolution);
        const std::vector<double> one = model.logLikelihood(9.0, 0.5);
        const st::LogLikelihoodMatrix many = model.logLikelihoods(
            std::vector<double>(kEvents, 9.0), std::vector<double>(kEvents, 0.5));
        const std::vector<double> combined = st::combineLogLikelihoods(many, log_prior_En);
        for (std::size_t i = 0; i < En.size(); ++i) {
            const double expected = log_prior_En[i] + (static_cast<double>(kEvents) * one[i]);
            CHECK_THAT(combined[i], WithinRel(expected, 1e-12) || WithinAbs(expected, 1e-9));
        }
    }

    SECTION("the prior is counted once") {
        const st::CaseAModel model = physicalModel(En, kToyResolution);
        const std::vector<double> single =
            st::singleEventLogPosterior(model, 9.0, 0.5, log_prior_En);
        const std::vector<double> combined = st::combineLogLikelihoods(
            model.logLikelihoods(std::vector<double>{9.0}, std::vector<double>{0.5}), log_prior_En);
        for (std::size_t i = 0; i < En.size(); ++i) {
            CHECK_THAT(combined[i], WithinULP(single[i], 0));
        }
    }
}

TEST_CASE("the sample overload reads the reco columns", "[stats][case_a]") {
    const st::Grid En = st::energyGrid(st::EnergyDomain{});
    const st::CaseAModel model = physicalModel(En, kToyResolution);
    riptide::EventSample sample;
    sample.reco_Ep = {4.0, 11.0};
    sample.reco_theta_p = {0.9, 0.2};
    const st::LogLikelihoodMatrix matrix = model.logLikelihoods(sample);
    REQUIRE(matrix.events == 2);
    const std::vector<double> second = model.logLikelihood(11.0, 0.2);
    for (std::size_t i = 0; i < En.size(); ++i) {
        CHECK_THAT(matrix.row(1)[i], WithinULP(second[i], 0));
    }
}

TEST_CASE("invalid Case A input is rejected", "[stats][case_a]") {
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    const st::Grid En = st::energyGrid(st::EnergyDomain{});
    const st::Grid theta = st::recoilAngleGrid();
    const std::vector<double> prior = st::physicalRecoilAngleLogPrior(theta);

    CHECK_THROWS_AS(st::CaseAModel(En, theta, prior, st::CaseAResolution{0.0, 0.08}),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::CaseAModel(En, theta, prior, st::CaseAResolution{0.10, -1.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::CaseAModel(En, theta, prior, st::CaseAResolution{kNaN, 0.08}),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::CaseAModel(En, theta, std::vector<double>(3, 0.0), kToyResolution),
                    std::invalid_argument);
    CHECK_THROWS_AS(
        st::CaseAModel(En, st::Grid(0.0, 2.0, 50), std::vector<double>(50, 0.0), kToyResolution),
        std::invalid_argument);

    const st::CaseAModel model(En, theta, prior, kToyResolution);
    CHECK_THROWS_AS(model.logLikelihood(kNaN, 0.5), std::invalid_argument);
    CHECK_THROWS_AS(model.logLikelihood(5.0, std::numeric_limits<double>::infinity()),
                    std::invalid_argument);
    CHECK_THROWS_AS(model.logLikelihoods(std::vector<double>{1.0}, std::vector<double>{}),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::singleEventLogPosterior(model, 5.0, 0.5, std::vector<double>(3, 0.0)),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::combineLogLikelihoods(
                        model.logLikelihoods(std::vector<double>{5.0}, std::vector<double>{0.5}),
                        std::vector<double>(3, 0.0)),
                    std::invalid_argument);
}

namespace {

// The likelihood summed over every theta cell, with no window: the reference
// the windowed row sum must reproduce.
std::vector<double> bruteForceLogLikelihood(const st::Grid& En, const st::Grid& theta,
                                            const std::vector<double>& log_prior_theta,
                                            const st::CaseAResolution& res, double Ep_hat,
                                            double theta_hat) {
    std::vector<double> result(En.size());
    std::vector<double> terms(theta.size());
    for (std::size_t i = 0; i < En.size(); ++i) {
        for (std::size_t j = 0; j < theta.size(); ++j) {
            const double zE = (Ep_hat - kin::protonEnergy(En[i], theta[j])) / res.sigma_Ep;
            const double zt = (theta_hat - theta[j]) / res.sigma_theta;
            terms[j] = (-0.5 * zE * zE) - (0.5 * zt * zt) + log_prior_theta[j];
        }
        result[i] = st::logSumExp(terms);
    }
    return result;
}

} // namespace

TEST_CASE("the theta window reproduces the sum over every cell", "[stats][case_a]") {
    const st::Grid En = st::energyGrid(st::EnergyDomain{});
    const st::Grid theta = st::recoilAngleGrid();
    // Observables across the whole range, including the edges the window
    // must handle: negative Ep_hat, Ep_hat above every Ep_max, theta_hat
    // outside [0, pi/2] and an event far from every cell.
    const std::vector<std::pair<double, double>> events = {{7.3, 0.55}, {0.02, 1.50}, {-0.3, 1.60},
                                                           {19.9, 0.0}, {25.0, 0.01}, {12.0, -0.2},
                                                           {3.0, 1.2},  {18.0, 1.4}};
    const std::vector<std::vector<double>> priors = {st::physicalRecoilAngleLogPrior(theta),
                                                     st::flatRecoilAngleLogPrior(theta)};
    const std::vector<st::CaseAResolution> resolutions = {kToyResolution, {0.02, 0.01}, {2.0, 0.5}};
    for (const auto& prior : priors) {
        for (const auto& res : resolutions) {
            const st::CaseAModel model(En, theta, prior, res);
            for (const auto& [Ep_hat, theta_hat] : events) {
                const std::vector<double> fast = model.logLikelihood(Ep_hat, theta_hat);
                const std::vector<double> ref =
                    bruteForceLogLikelihood(En, theta, prior, res, Ep_hat, theta_hat);
                double worst = 0.0;
                for (std::size_t i = 0; i < ref.size(); ++i) {
                    worst = std::max(worst, std::fabs(fast[i] - ref[i]) /
                                                std::max(1.0, std::fabs(ref[i])));
                }
                // The skipped mass is below 1e-19; what is left is the
                // rounding of the kept terms.
                CHECK(worst < 1e-13);
            }
        }
    }
}
