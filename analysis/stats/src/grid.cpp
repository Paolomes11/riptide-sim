#include "riptide/stats/grid.hpp"

#include <cmath>
#include <stdexcept>

namespace riptide::stats {

namespace {
constexpr std::size_t kMinGridPoints = 2;
} // namespace

Grid::Grid(double lo, double hi, std::size_t n) {
    if (!std::isfinite(lo) || !std::isfinite(hi) || !(lo < hi)) {
        throw std::invalid_argument("grid edges must be finite with lo < hi");
    }
    if (n < kMinGridPoints) {
        throw std::invalid_argument("a grid needs at least 2 points");
    }
    step_ = (hi - lo) / static_cast<double>(n - 1);
    values_.resize(n);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        values_[i] = lo + static_cast<double>(i) * step_;
    }
    values_.back() = hi;
}

Grid energyGrid(const EnergyDomain& domain) {
    return {domain.En_min, domain.En_max, domain.points};
}

EnergyDomain sameStepDomain(const EnergyDomain& base, double En_max) {
    const Grid reference = energyGrid(base);
    if (!std::isfinite(En_max) || !(En_max > base.En_min)) {
        throw std::invalid_argument("the widened domain must end above its lower edge");
    }
    const double intervals = std::round((En_max - base.En_min) / reference.step());
    EnergyDomain widened = base;
    widened.En_max = En_max;
    widened.points = static_cast<std::size_t>(intervals) + 1;
    return widened;
}

Grid recoilAngleGrid(std::size_t points) {
    return {0.0, 0.5 * units::pi, points};
}

} // namespace riptide::stats
