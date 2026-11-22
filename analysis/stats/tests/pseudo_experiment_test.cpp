#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment.hpp"
#include "riptide/stats/validate.hpp"

namespace st = riptide::stats;
namespace kin = riptide::kinematics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;

namespace {

constexpr std::uint64_t kSeed = 20260907;
constexpr st::CaseAResolution kToyResolution{0.10, 0.08};
constexpr double kMcSigmas = 4.0;

struct Moments {
    double mean = 0.0;
    double variance = 0.0;
};

Moments moments(const std::vector<double>& x) {
    double sum = 0.0;
    for (const double v : x) {
        sum += v;
    }
    const double mean = sum / static_cast<double>(x.size());
    double squares = 0.0;
    for (const double v : x) {
        squares += (v - mean) * (v - mean);
    }
    return {mean, squares / static_cast<double>(x.size())};
}

} // namespace

TEST_CASE("pseudo-events have Ep uniform below Ep_max", "[stats][pseudo]") {
    constexpr std::size_t kSamples = 1000000;
    constexpr double kEn = 14.0;
    const auto sample = st::pseudoSampleAtEnergy(kEn, kSamples, kToyResolution, kSeed);
    const double width = kin::maxProtonEnergy(kEn);
    const auto m = moments(sample.truth_Ep);
    const double n = static_cast<double>(kSamples);
    // Uniform on [0, L): mean L/2, variance L^2/12, and the variance of
    // (x - mu)^2 is L^4/80 - L^4/144 = L^4/180.
    CHECK_THAT(m.mean, WithinAbs(0.5 * width, kMcSigmas * width / std::sqrt(12.0 * n)));
    CHECK_THAT(m.variance,
               WithinAbs(width * width / 12.0, kMcSigmas * width * width / std::sqrt(180.0 * n)));
    std::size_t outside = 0;
    for (std::size_t k = 0; k < kSamples; ++k) {
        const bool in_range = sample.truth_Ep[k] >= 0.0 && sample.truth_Ep[k] < width &&
                              std::fabs(sample.truth_En[k] - kEn) <= 0.0;
        outside += in_range ? 0U : 1U;
    }
    CHECK(outside == 0);
}

TEST_CASE("the recoil angle has density sin(2 theta) at low energy", "[stats][pseudo]") {
    // Non-relativistically cos^2 theta_p = Ep / En with Ep uniform, so
    // sin^2 theta_p is uniform on [0, 1): mean 1/2, variance 1/12. The
    // relativistic and mass corrections are ~En / m, negligible at 10 keV.
    constexpr std::size_t kSamples = 1000000;
    constexpr double kEn = 0.01;
    const auto sample = st::pseudoSampleAtEnergy(kEn, kSamples, kToyResolution, kSeed);
    constexpr double kHalfPi = 0.5 * riptide::units::pi;
    std::vector<double> sin2(kSamples);
    std::size_t outside = 0;
    for (std::size_t k = 0; k < kSamples; ++k) {
        const double s = std::sin(sample.truth_theta_p[k]);
        sin2[k] = s * s;
        outside += (sample.truth_theta_p[k] >= 0.0 && sample.truth_theta_p[k] <= kHalfPi) ? 0U : 1U;
    }
    CHECK(outside == 0);
    const auto m = moments(sin2);
    const double n = static_cast<double>(kSamples);
    CHECK_THAT(m.mean, WithinAbs(0.5, kMcSigmas / std::sqrt(12.0 * n)));
    CHECK_THAT(m.variance, WithinAbs(1.0 / 12.0, kMcSigmas / std::sqrt(180.0 * n)));
}

TEST_CASE("truth Ep and theta_p follow the exact kinematics", "[stats][pseudo]") {
    const auto sample = st::pseudoSampleAtEnergy(17.0, 1000, kToyResolution, kSeed);
    for (std::size_t k = 0; k < sample.size(); ++k) {
        CHECK_THAT(kin::protonEnergy(17.0, sample.truth_theta_p[k]),
                   WithinAbs(sample.truth_Ep[k], 1e-9));
    }
}

TEST_CASE("measurement noise has the declared widths", "[stats][pseudo]") {
    constexpr std::size_t kSamples = 200000;
    const auto sample = st::pseudoSampleAtEnergy(14.0, kSamples, kToyResolution, kSeed);
    std::vector<double> dEp(kSamples);
    std::vector<double> dtheta(kSamples);
    for (std::size_t k = 0; k < kSamples; ++k) {
        dEp[k] = sample.reco_Ep[k] - sample.truth_Ep[k];
        dtheta[k] = sample.reco_theta_p[k] - sample.truth_theta_p[k];
    }
    const double n = static_cast<double>(kSamples);
    const auto mEp = moments(dEp);
    const auto mth = moments(dtheta);
    const double sEp = kToyResolution.sigma_Ep;
    const double sth = kToyResolution.sigma_theta;
    CHECK_THAT(mEp.mean, WithinAbs(0.0, kMcSigmas * sEp / std::sqrt(n)));
    CHECK_THAT(std::sqrt(mEp.variance), WithinAbs(sEp, kMcSigmas * sEp / std::sqrt(2.0 * n)));
    CHECK_THAT(mth.mean, WithinAbs(0.0, kMcSigmas * sth / std::sqrt(n)));
    CHECK_THAT(std::sqrt(mth.variance), WithinAbs(sth, kMcSigmas * sth / std::sqrt(2.0 * n)));
}

TEST_CASE("pseudo-samples are reproducible event by event", "[stats][pseudo]") {
    const st::EnergyDomain domain;
    const auto a = st::pseudoSampleFromPrior(domain, 500, kToyResolution, kSeed);
    const auto b = st::pseudoSampleFromPrior(domain, 500, kToyResolution, kSeed);
    // The tail of a sample is the sample that starts at its first event id.
    const auto tail = st::pseudoSampleFromPrior(domain, 200, kToyResolution, kSeed, 300);
    const auto other = st::pseudoSampleFromPrior(domain, 500, kToyResolution, kSeed + 1);
    std::size_t differing = 0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        REQUIRE_THAT(a.truth_En[k], WithinULP(b.truth_En[k], 0));
        REQUIRE_THAT(a.reco_Ep[k], WithinULP(b.reco_Ep[k], 0));
        REQUIRE_THAT(a.reco_theta_p[k], WithinULP(b.reco_theta_p[k], 0));
        REQUIRE(a.truth_En[k] >= domain.En_min);
        REQUIRE(a.truth_En[k] < domain.En_max);
        differing += std::fabs(a.reco_Ep[k] - other.reco_Ep[k]) > 0.0 ? 1U : 0U;
    }
    for (std::size_t k = 0; k < tail.size(); ++k) {
        REQUIRE_THAT(tail.truth_En[k], WithinULP(a.truth_En[300 + k], 0));
        REQUIRE_THAT(tail.reco_Ep[k], WithinULP(a.reco_Ep[300 + k], 0));
        REQUIRE_THAT(tail.reco_theta_p[k], WithinULP(a.reco_theta_p[300 + k], 0));
    }
    CHECK(differing == a.size());
}

TEST_CASE("pseudo-experiment inputs are checked", "[stats][pseudo]") {
    CHECK_THROWS_AS(st::pseudoSampleAtEnergy(14.0, 1, {0.0, 0.08}, kSeed), std::invalid_argument);
    CHECK_THROWS_AS(st::pseudoSampleAtEnergy(-1.0, 1, kToyResolution, kSeed),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::pseudoSampleFromPrior({20.0, 10.0, 500}, 1, kToyResolution, kSeed),
                    std::invalid_argument);
}

TEST_CASE("posterior A is calibrated when the truth comes from the prior",
          "[stats][pseudo][calibration]") {
    // 2000 events: the binomial error at 68% is sqrt(0.68 * 0.32 / 2000) =
    // 0.010; the threshold is 3 errors, as in the toy.
    constexpr std::size_t kEvents = 2000;
    const std::vector<double> levels{0.68, 0.90};
    const st::EnergyDomain domain;
    const st::Grid En = st::energyGrid(domain);
    const st::Grid theta = st::recoilAngleGrid();
    const st::CaseAModel model(En, theta, st::physicalRecoilAngleLogPrior(theta), kToyResolution);
    const std::vector<double> log_prior = st::flatEnergyLogPrior(En);

    const auto sample = st::pseudoSampleFromPrior(domain, kEvents, kToyResolution, kSeed);
    const st::LogLikelihoodMatrix loglik =
        model.logLikelihoods(sample.reco_Ep, sample.reco_theta_p);
    st::IntervalTable table;
    table.levels = levels;
    table.values.reserve(kEvents * levels.size());
    std::vector<double> log_post(En.size());
    for (std::size_t k = 0; k < kEvents; ++k) {
        const double* row = loglik.row(k);
        for (std::size_t i = 0; i < En.size(); ++i) {
            log_post[i] = row[i] + log_prior[i];
        }
        for (const auto& interval : st::credibleIntervals(En, log_post.data(), levels)) {
            table.values.push_back(interval);
        }
    }
    const std::vector<double> c = st::coverage(sample.truth_En, table);
    for (std::size_t l = 0; l < levels.size(); ++l) {
        const double binomial =
            std::sqrt(levels[l] * (1.0 - levels[l]) / static_cast<double>(kEvents));
        CHECK_THAT(c[l], WithinAbs(levels[l], 3.0 * binomial));
    }
}
