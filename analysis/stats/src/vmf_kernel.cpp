#include "riptide/stats/vmf_kernel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "riptide/common/units.hpp"

namespace riptide::stats {
namespace {

// The asymptotic series is stopped once a term is below this fraction of the
// sum: well under one ulp. Above kBesselAsymptoticFrom the terms keep falling
// for k < 2x, far past where this is reached (k ~ 10 at x = 30).
constexpr double kSeriesRelativeTolerance = 1e-17;
constexpr int kSeriesMaxTerms = 40;

double logI0eAsymptotic(double x) {
    const double eight_x = 8.0 * x;
    double term = 1.0;
    double sum = 1.0;
    for (int k = 1; k <= kSeriesMaxTerms; ++k) {
        const auto order = static_cast<double>(k);
        const double odd = (2.0 * order) - 1.0;
        term *= odd * odd / (order * eight_x);
        sum += term;
        if (term < kSeriesRelativeTolerance * sum) {
            break;
        }
    }
    return std::log(sum) - (0.5 * std::log(2.0 * units::pi * x));
}

} // namespace

double logI0e(double x) {
    if (std::isnan(x) || x < 0.0) {
        throw std::invalid_argument("logI0e: the argument must be >= 0");
    }
    if (std::isinf(x)) {
        return -std::numeric_limits<double>::infinity();
    }
    if (x > kBesselAsymptoticFrom) {
        return logI0eAsymptotic(x);
    }
    // I0(x) e^{-x} at x <= 30 is ~1e13 at most before the product: no
    // overflow, and the exponential is applied in the log.
    return std::log(std::cyl_bessel_i(0.0, x)) - x;
}

double vmfConcentration(double sigma_theta) {
    if (!std::isfinite(sigma_theta) || sigma_theta <= 0.0) {
        throw std::invalid_argument("vmfConcentration: sigma_theta must be finite and > 0");
    }
    return 1.0 / (sigma_theta * sigma_theta);
}

double logMeridianVmfKernel(double theta_obs, double theta, double kappa) noexcept {
    const double sin_theta = std::sin(theta);
    // cos(d) - 1 = -2 sin^2(d/2) has no cancellation for d << 1, where a large
    // kappa puts all the weight.
    const double half_difference = std::sin(0.5 * (theta_obs - theta));
    const double log_norm = std::log(kappa) - std::log1p(-std::exp(-2.0 * kappa));
    // The argument of I0e is >= 0 on [0, pi]; the max only absorbs an angle
    // rounded just past pi, whose sine is a tiny negative number.
    const double bessel_argument = std::max(0.0, kappa * std::sin(theta_obs) * sin_theta);
    return log_norm + std::log(sin_theta) - (2.0 * kappa * half_difference * half_difference) +
           logI0e(bessel_argument);
}

double logTrackDensity(double theta) noexcept {
    return std::log(std::max(0.0, std::cos(theta)) / units::pi);
}

} // namespace riptide::stats
