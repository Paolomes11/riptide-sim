#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

#include "riptide/common/units.hpp"
#include "riptide/stats/vmf_kernel.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::stats::logI0e;
using riptide::stats::logMeridianVmfKernel;
using riptide::stats::logTrackDensity;
using riptide::stats::vmfConcentration;
using riptide::units::pi;

struct BesselReference {
    double x;
    double i0e;
};

// scipy.special.i0e (scipy 1.18.1), printed with repr: exp(-x) I0(x).
constexpr BesselReference kI0eReference[] = {
    {0.0, 1.0},
    {1e-08, 0.9999999900000003},
    {0.001, 0.9990007495835156},
    {0.1, 0.9071009257823011},
    {0.5, 0.64503527044915},
    {1.0, 0.46575960759364043},
    {2.5, 0.27004644161220276},
    {7.0, 0.15373774467288123},
    {15.0, 0.1038995314488227},
    {30.0, 0.0731459464822373},
    {49.0, 0.057138847093590205},
    {50.0, 0.056561626647454184},
    {51.0, 0.05600155428684504},
    {80.0, 0.04467329178227527},
    {156.25, 0.03194100727322136},
    {500.0, 0.017845706500153165},
    {10000.0, 0.0039894726746047314},
    {1000000.0, 0.00039894233026924577},
};

// int_0^pi sin(theta_obs) K(theta_obs | theta) dtheta_obs by trapezoids.
double meridianIntegral(double theta, double kappa, std::size_t points) {
    const double step = pi / static_cast<double>(points - 1);
    double sum = 0.0;
    for (std::size_t i = 0; i < points; ++i) {
        const double theta_obs = step * static_cast<double>(i);
        const double weight = (i == 0 || i + 1 == points) ? 0.5 : 1.0;
        sum +=
            weight * std::sin(theta_obs) * std::exp(logMeridianVmfKernel(theta_obs, theta, kappa));
    }
    return sum * step;
}

} // namespace

TEST_CASE("logI0e matches the reference values", "[vmf_kernel]") {
    for (const auto& reference : kI0eReference) {
        CAPTURE(reference.x);
        CHECK_THAT(std::exp(logI0e(reference.x)), WithinRel(reference.i0e, 1e-12));
    }
}

TEST_CASE("logI0e is continuous at the switch to the asymptotic series", "[vmf_kernel]") {
    const double below = std::nextafter(riptide::stats::kBesselAsymptoticFrom, 0.0);
    const double above = std::nextafter(riptide::stats::kBesselAsymptoticFrom, 100.0);
    CHECK_THAT(logI0e(above), WithinAbs(logI0e(below), 1e-14));
}

TEST_CASE("logI0e edge cases", "[vmf_kernel]") {
    CHECK(logI0e(std::numeric_limits<double>::infinity()) < 0.0);
    CHECK(std::isinf(logI0e(std::numeric_limits<double>::infinity())));
    CHECK_THROWS_AS(logI0e(-1e-300), std::invalid_argument);
    CHECK_THROWS_AS(logI0e(std::nan("")), std::invalid_argument);
}

TEST_CASE("vmfConcentration is 1 / sigma^2", "[vmf_kernel]") {
    CHECK_THAT(vmfConcentration(0.08), WithinRel(156.25, 1e-15));
    CHECK_THROWS_AS(vmfConcentration(0.0), std::invalid_argument);
    CHECK_THROWS_AS(vmfConcentration(-0.1), std::invalid_argument);
    CHECK_THROWS_AS(vmfConcentration(std::numeric_limits<double>::infinity()),
                    std::invalid_argument);
}

TEST_CASE("the meridian kernel is normalised over the observed track", "[vmf_kernel]") {
    // The vMF density integrates to 1 over the observed direction, and the
    // azimuth of the true track contributes its circle 2 pi sin(theta): the
    // integral over theta_obs is sin(theta), also next to the pole, where
    // part of the density lies past it.
    constexpr std::size_t kPoints = 40001;
    for (const double sigma : {0.3, 0.08, 0.02}) {
        const double kappa = vmfConcentration(sigma);
        for (const double theta : {0.03, 0.3, 0.8, 1.3, riptide::units::pi / 2.0}) {
            CAPTURE(sigma, theta);
            CHECK_THAT(meridianIntegral(theta, kappa, kPoints) / std::sin(theta),
                       WithinRel(1.0, 1e-6));
        }
    }
}

TEST_CASE("the meridian kernel has the Gaussian limit at large kappa", "[vmf_kernel]") {
    // I0e(x) -> 1 / sqrt(2 pi x): K -> sqrt(sin theta / sin theta_obs)
    // N(theta_obs; theta, sigma). The residual is the 1/(8x) term (~3e-7 here)
    // and the quartic term of 2 sin^2(d/2) (~4e-6 at 3 sigma).
    constexpr double kSigma = 1e-3;
    constexpr double kTheta = 0.7;
    const double kappa = vmfConcentration(kSigma);
    for (const double pulls : {-3.0, -1.0, 0.0, 0.5, 2.0, 3.0}) {
        const double theta_obs = kTheta + (pulls * kSigma);
        const double log_gauss = (-0.5 * pulls * pulls) - std::log(kSigma * std::sqrt(2.0 * pi));
        const double log_limit =
            log_gauss + (0.5 * std::log(std::sin(kTheta) / std::sin(theta_obs)));
        CAPTURE(pulls);
        CHECK_THAT(logMeridianVmfKernel(theta_obs, kTheta, kappa), WithinAbs(log_limit, 1e-5));
    }
}

TEST_CASE("the meridian kernel vanishes at theta = 0", "[vmf_kernel]") {
    const double value = logMeridianVmfKernel(0.1, 0.0, vmfConcentration(0.08));
    CHECK(std::isinf(value));
    CHECK(value < 0.0);
    CHECK(std::isfinite(logMeridianVmfKernel(pi, 0.5, vmfConcentration(0.08))));
}

TEST_CASE("logTrackDensity is log(cos(theta) / pi) on the forward hemisphere", "[vmf_kernel]") {
    CHECK_THAT(logTrackDensity(0.0), WithinAbs(-std::log(pi), 1e-15));
    CHECK_THAT(logTrackDensity(1.0), WithinAbs(std::log(std::cos(1.0) / pi), 1e-15));
    CHECK(logTrackDensity(pi / 2.0) < -35.0);
    const double backward = logTrackDensity(2.0);
    CHECK(std::isinf(backward));
    CHECK(backward < 0.0);
}
