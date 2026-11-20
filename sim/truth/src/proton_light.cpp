#include "riptide/sim/proton_light.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <string>
#include <utility>

namespace riptide::sim {

namespace {

// Power law through (x0, y0) and (x1, y1), all positive, evaluated at x.
double logLog(double x0, double y0, double x1, double y1, double x) {
    const double slope = std::log(y1 / y0) / std::log(x1 / x0);
    return y0 * std::exp(slope * std::log(x / x0));
}

bool finitePositive(double v) {
    return std::isfinite(v) && v > 0.0;
}

} // namespace

ProtonLightTable::ProtonLightTable(std::vector<double> K, std::vector<double> S, double kB)
    : K_(std::move(K)), S_(std::move(S)), kB_(kB) {
    if (K_.size() < 2 || K_.size() != S_.size()) {
        throw ProtonLightError(
            "proton light: need at least two energies, each with a stopping power");
    }
    if (!(kB_ >= 0.0) || !std::isfinite(kB_)) {
        throw ProtonLightError("proton light: kB must be finite and not negative");
    }
    for (std::size_t i = 0; i < K_.size(); ++i) {
        if (!finitePositive(K_[i]) || !finitePositive(S_[i])) {
            throw ProtonLightError("proton light: energies and stopping powers must be finite "
                                   "and positive");
        }
        if (i > 0 && !(K_[i] > K_[i - 1])) {
            throw ProtonLightError("proton light: energies must increase strictly");
        }
    }
    cumulative_.reserve(K_.size());
    cumulative_.push_back(K_.front() * quench(S_.front()));
    for (std::size_t i = 0; i + 1 < K_.size(); ++i) {
        cumulative_.push_back(cumulative_.back() + intervalLight(i, K_[i + 1]));
    }
}

double ProtonLightTable::quench(double S) const {
    return 1.0 / (1.0 + kB_ * S);
}

// Simpson's rule from K[i] to upper <= K[i+1], with S log-log in between, so
// the integrand is smooth inside every panel.
double ProtonLightTable::intervalLight(std::size_t i, double upper) const {
    const double e0 = K_[i];
    const double e1 = K_[i + 1];
    const double s0 = S_[i];
    const double s1 = S_[i + 1];
    const auto integrand = [&](double E) { return quench(logLog(e0, s0, e1, s1, E)); };
    const double width = (upper - e0) / kProtonLightPanels;
    double sum = 0.0;
    for (int k = 0; k < kProtonLightPanels; ++k) {
        const double a = e0 + k * width;
        const double b = (k + 1 == kProtonLightPanels) ? upper : a + width;
        sum += (b - a) * (integrand(a) + 4.0 * integrand(0.5 * (a + b)) + integrand(b));
    }
    return sum / 6.0;
}

double ProtonLightTable::light(double K) const {
    if (!(K >= 0.0) || !(K <= K_.back())) {
        throw ProtonLightError("proton light: energy " + std::to_string(K) +
                               " MeV negative or above the table");
    }
    if (!(K > K_.front())) {
        return K * quench(S_.front());
    }
    // K in (K[lo], K[lo+1]]: the full intervals below lo are in the sum.
    const auto upper = std::lower_bound(K_.begin(), K_.end(), K);
    const auto lo = static_cast<std::size_t>(std::distance(K_.begin(), upper)) - 1;
    return cumulative_[lo] + intervalLight(lo, K);
}

double ProtonLightTable::visibleEnergy(double K_pre, double edep) const {
    if (!(edep >= 0.0) || !std::isfinite(edep)) {
        throw ProtonLightError("proton light: deposit must be finite and not negative");
    }
    const double continuous = std::min(edep, K_pre);
    const double excess = edep - continuous;
    const double visible = light(K_pre) - light(K_pre - continuous) + excess * quench(S_.front());
    return std::clamp(visible, 0.0, edep);
}

} // namespace riptide::sim
