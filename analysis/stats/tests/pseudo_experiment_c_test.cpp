#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

#include "riptide/common/direction.hpp"
#include "riptide/common/random.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_c.hpp"

namespace st = riptide::stats;
namespace geo = riptide::geometry;
namespace rnd = riptide::random;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace {

constexpr std::uint64_t kSeed = 20260907;
const st::EnergyDomain kDomain{10.0, 20.0, st::kEnergyGridPoints};
const st::CaseBResolution kResolution{0.10, 0.08};

// Bound of a binomial fraction at 3 standard deviations.
double threeSigma(double p, std::size_t n) {
    return 3.0 * std::sqrt(p * (1.0 - p) / static_cast<double>(n));
}

double standardNormalDensity(double z) {
    return std::exp(-0.5 * z * z) / std::sqrt(2.0 * 3.14159265358979323846);
}

} // namespace

TEST_CASE("the mass of N(mu, sigma) within one sigma is erf(1 / sqrt 2)", "[case_c_pseudo]") {
    const st::EnergyDomain one_sigma{14.0, 16.0, st::kEnergyGridPoints};
    CHECK_THAT(st::normalMassOnDomain(one_sigma, 15.0, 1.0),
               WithinRel(std::erf(1.0 / std::sqrt(2.0)), 1e-14));
    // Far tails on both sides: Q(10) - Q(20) = erfc(10 / sqrt 2) / 2 to
    // double precision, the same mirrored.
    const double tail = 0.5 * std::erfc(10.0 / std::sqrt(2.0));
    CHECK_THAT(st::normalMassOnDomain(kDomain, 0.0, 1.0), WithinRel(tail, 1e-14));
    const st::EnergyDomain below{1.0, 11.0, st::kEnergyGridPoints};
    CHECK_THAT(st::normalMassOnDomain(below, 21.0, 1.0), WithinRel(tail, 1e-14));
}

TEST_CASE("Case C truths stay in their ranges", "[case_c_pseudo]") {
    for (std::uint64_t i = 0; i < 200; ++i) {
        auto engine = rnd::perEventEngine(kSeed, i);
        const st::CaseCTruth truth = st::drawCaseCTruth(st::kCaseCChecklistTruth, engine);
        CHECK(truth.mu >= 13.0);
        CHECK(truth.mu <= 17.0);
        CHECK(truth.sigma >= 0.8 * (1.0 - 1e-15));
        CHECK(truth.sigma <= 2.4 * (1.0 + 1e-15));
        const geo::Direction& w = truth.omega_n;
        CHECK_THAT(std::sqrt((w.x * w.x) + (w.y * w.y) + (w.z * w.z)), WithinAbs(1.0, 1e-15));
    }
    auto engine = rnd::perEventEngine(kSeed, 0);
    CHECK_THROWS_AS(st::drawCaseCTruth({13.0, 17.0, 0.0, 2.4}, engine), std::invalid_argument);
    CHECK_THROWS_AS(st::drawCaseCTruth({17.0, 13.0, 0.8, 2.4}, engine), std::invalid_argument);
}

TEST_CASE("truncated edge experiments put no event outside the domain", "[case_c_pseudo]") {
    std::size_t events = 0;
    for (std::uint64_t i = 0; i < 20; ++i) {
        const st::CaseCPseudoExperiment experiment =
            st::caseCPseudoExperiment(kDomain, 200, kResolution, st::kCaseCEdgeTruth,
                                      st::CaseCSpectrum::truncated, kSeed, 1, i);
        CHECK(experiment.truth.mu >= 17.5);
        CHECK(experiment.truth.mu <= 19.0);
        for (const double En : experiment.sample.truth_En) {
            CHECK(En >= kDomain.En_min);
            CHECK(En <= kDomain.En_max);
        }
        events += experiment.sample.size();
    }
    CHECK(events == 4000);
}

TEST_CASE("truncated energies have the mean of the truncated Gaussian", "[case_c_pseudo]") {
    // mu = 19, sigma = 2: the domain cuts the spectrum at +0.5 sigma.
    constexpr double kMu = 19.0;
    constexpr double kSigma = 2.0;
    constexpr std::size_t kDraws = 20000;
    double sum = 0.0;
    double sum_sq = 0.0;
    for (std::size_t k = 0; k < kDraws; ++k) {
        auto engine = rnd::perEventEngine(kSeed, k);
        const double En = st::drawTruncatedNormalEnergy(kDomain, kMu, kSigma, engine);
        sum += En;
        sum_sq += En * En;
    }
    const double a = (kDomain.En_min - kMu) / kSigma;
    const double b = (kDomain.En_max - kMu) / kSigma;
    const double mass = st::normalMassOnDomain(kDomain, kMu, kSigma);
    const double expected =
        kMu + (kSigma * (standardNormalDensity(a) - standardNormalDensity(b)) / mass);
    const double mean = sum / static_cast<double>(kDraws);
    const double sd = std::sqrt((sum_sq / static_cast<double>(kDraws)) - (mean * mean));
    CHECK(std::abs(mean - expected) < 4.0 * sd / std::sqrt(static_cast<double>(kDraws)));
}

TEST_CASE("untruncated energies fall outside the domain at the Gaussian rate", "[case_c_pseudo]") {
    // mu = 13, sigma = 2.4 on [10, 20]: about 10.7% outside, mostly below.
    const st::CaseCTruth truth{13.0, 2.4, geo::Direction{0.0, 0.0, 1.0}};
    constexpr std::size_t kDraws = 20000;
    std::size_t outside = 0;
    for (std::size_t k = 0; k < kDraws; ++k) {
        auto engine = rnd::perEventEngine(kSeed, k);
        const st::CaseBPseudoEvent event = st::drawCaseCPseudoEvent(
            kDomain, truth, st::CaseCSpectrum::untruncated, kResolution, engine);
        if (event.truth_En < kDomain.En_min || event.truth_En > kDomain.En_max) {
            ++outside;
        }
    }
    const double p_out = 1.0 - st::normalMassOnDomain(kDomain, truth.mu, truth.sigma);
    CHECK_THAT(p_out, WithinAbs(0.1074, 1e-3));
    const double fraction = static_cast<double>(outside) / static_cast<double>(kDraws);
    CHECK(std::abs(fraction - p_out) < threeSigma(p_out, kDraws));
}

TEST_CASE("Case C events do not depend on how many are drawn", "[case_c_pseudo]") {
    const st::CaseCPseudoExperiment small =
        st::caseCPseudoExperiment(kDomain, 10, kResolution, st::kCaseCChecklistTruth,
                                  st::CaseCSpectrum::truncated, kSeed, 2, 7);
    const st::CaseCPseudoExperiment large =
        st::caseCPseudoExperiment(kDomain, 30, kResolution, st::kCaseCChecklistTruth,
                                  st::CaseCSpectrum::truncated, kSeed, 2, 7);
    CHECK_THAT(small.truth.mu, WithinULP(large.truth.mu, 0));
    CHECK_THAT(small.truth.sigma, WithinULP(large.truth.sigma, 0));
    CHECK_THAT(small.truth.omega_n.x, WithinULP(large.truth.omega_n.x, 0));
    CHECK_THAT(small.truth.omega_n.z, WithinULP(large.truth.omega_n.z, 0));
    for (std::size_t k = 0; k < small.sample.size(); ++k) {
        CHECK_THAT(small.sample.truth_En[k], WithinULP(large.sample.truth_En[k], 0));
        CHECK_THAT(small.sample.reco_Ep[k], WithinULP(large.sample.reco_Ep[k], 0));
        CHECK_THAT(small.sample.reco_track[k].y, WithinULP(large.sample.reco_track[k].y, 0));
    }
}

TEST_CASE("Case C generation rejects spectra it cannot draw", "[case_c_pseudo]") {
    auto engine = rnd::perEventEngine(kSeed, 0);
    CHECK_THROWS_AS(st::drawTruncatedNormalEnergy(kDomain, 40.0, 1.0, engine),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::drawTruncatedNormalEnergy(kDomain, 15.0, 0.0, engine),
                    std::invalid_argument);
    const st::CaseCTruth negative{-5.0, 0.1, geo::Direction{0.0, 0.0, 1.0}};
    CHECK_THROWS_AS(st::drawCaseCPseudoEvent(kDomain, negative, st::CaseCSpectrum::untruncated,
                                             kResolution, engine),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::caseCPseudoExperiment(kDomain, 10, {0.0, 0.08}, st::kCaseCChecklistTruth,
                                              st::CaseCSpectrum::truncated, kSeed, 0, 0),
                    std::invalid_argument);
}
