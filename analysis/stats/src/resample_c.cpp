#include "riptide/stats/resample_c.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/stats/split.hpp"

namespace riptide::stats {

namespace {

// Relative slack allowed when checking that the bins tile the domain: the
// width over the step is an integer up to the rounding of the two numbers.
constexpr double kTilingTolerance = 1e-9;

std::size_t binCount(const EnergyDomain& domain, double step) {
    if (!(step > 0.0) || !std::isfinite(step)) {
        throw std::invalid_argument("the bin width must be positive and finite");
    }
    if (!(domain.En_max > domain.En_min) || !std::isfinite(domain.En_min) ||
        !std::isfinite(domain.En_max)) {
        throw std::invalid_argument("the density needs a finite domain with En_min < En_max");
    }
    const double ratio = (domain.En_max - domain.En_min) / step;
    const double rounded = std::round(ratio);
    if (rounded < 1.0 || std::abs(ratio - rounded) > kTilingTolerance * ratio) {
        throw std::invalid_argument("the bin width does not tile the domain");
    }
    return static_cast<std::size_t>(rounded);
}

// Unnormalised Gaussian kernel exp(-z^2 / 2): the constants cancel in f / f.
double gaussianKernel(double En, double mu, double sigma) {
    const double z = (En - mu) / sigma;
    return std::exp(-0.5 * z * z);
}

double normaliser(const EnergyDomain& domain, double mu, double sigma) {
    const double mass = normalMassOnDomain(domain, mu, sigma);
    if (mass < kMinTruncatedMass) {
        throw std::invalid_argument("the Case C spectrum has too little mass on the domain");
    }
    return 1.0 / (sigma * std::sqrt(2.0 * units::pi) * mass);
}

} // namespace

std::size_t BinnedDensity::binOf(double En) const {
    if (!(En >= domain.En_min && En <= domain.En_max)) {
        throw std::invalid_argument("energy outside the domain of the binned density");
    }
    const auto bin = static_cast<std::size_t>((En - domain.En_min) / step);
    return std::min(bin, density.size() - 1);
}

BinnedDensity histogramDensity(const std::vector<double>& values, const EnergyDomain& domain,
                               double step) {
    BinnedDensity out;
    out.domain = domain;
    out.step = step;
    out.density.assign(binCount(domain, step), 0.0);
    for (const double En : values) {
        out.density[out.binOf(En)] += 1.0;
    }
    const double scale = 1.0 / (static_cast<double>(values.size()) * step);
    for (double& d : out.density) {
        if (!(d > 0.0)) {
            throw std::invalid_argument("an empty bin in the pool: c would be infinite");
        }
        d *= scale;
    }
    return out;
}

double truncatedNormalDensity(const EnergyDomain& domain, double mu, double sigma, double En) {
    const double norm = normaliser(domain, mu, sigma);
    if (!(En >= domain.En_min && En <= domain.En_max)) {
        return 0.0;
    }
    return norm * gaussianKernel(En, mu, sigma);
}

double resampleBound(const BinnedDensity& density, double mu, double sigma) {
    const double norm = normaliser(density.domain, mu, sigma);
    double bound = 0.0;
    for (std::size_t b = 0; b < density.density.size(); ++b) {
        const double lo = density.domain.En_min + (static_cast<double>(b) * density.step);
        const double hi = std::min(lo + density.step, density.domain.En_max);
        const double peak = std::clamp(mu, lo, hi);
        bound = std::max(bound, norm * gaussianKernel(peak, mu, sigma) / density.density[b]);
    }
    return bound;
}

ResampledRow resampleRow(const std::vector<double>& truth_En, const BinnedDensity& density,
                         const std::vector<CaseCTruth>& truths, std::size_t events,
                         random::Engine& engine) {
    if (events == 0) {
        throw std::invalid_argument("a resampled experiment needs at least one event");
    }
    ResampledRow row;
    row.planned = truths.size();
    if (truth_En.empty() || truths.empty()) {
        return row;
    }
    const std::vector<std::size_t> order =
        disjointExperiments(truth_En.size(), truth_En.size(), 1, engine).front();
    std::size_t next = 0;
    for (const CaseCTruth& truth : truths) {
        const double bound = resampleBound(density, truth.mu, truth.sigma);
        const double norm = normaliser(density.domain, truth.mu, truth.sigma);
        std::vector<std::size_t> kept;
        kept.reserve(events);
        const std::size_t start = next;
        while (kept.size() < events && next < order.size()) {
            const std::size_t index = order[next++];
            const double En = truth_En[index];
            const double g = density.density[density.binOf(En)];
            const double f = norm * gaussianKernel(En, truth.mu, truth.sigma);
            // u < f / (c g) written without the division; f <= c g by the
            // choice of c, so the probability never exceeds 1.
            if (random::uniform53(engine) * bound * g < f) {
                kept.push_back(index);
            }
        }
        row.visited_total = next;
        if (kept.size() < events) {
            break;
        }
        row.experiments.push_back(std::move(kept));
        row.visited.push_back(next - start);
        row.bounds.push_back(bound);
    }
    return row;
}

} // namespace riptide::stats
