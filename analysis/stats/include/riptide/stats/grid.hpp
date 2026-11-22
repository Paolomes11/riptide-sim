#pragma once

// Uniform grids of the Case A posterior (D32): E_n on a parametric domain and
// the proton angle theta_p on [0, pi/2]. Grids are contiguous doubles fixed at
// construction, so every function reading them can assume the same layout.

#include <cstddef>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// Default grid sizes, the same as the toy model so its numbers can be compared
// at the same discretisation.
inline constexpr std::size_t kEnergyGridPoints = 500;
inline constexpr std::size_t kRecoilAngleGridPoints = 500;

class Grid {
public:
    // n points from lo to hi included, like numpy.linspace: point i is
    // lo + i * step and the last one is hi exactly. Throws
    // std::invalid_argument unless lo and hi are finite, lo < hi and n >= 2.
    Grid(double lo, double hi, std::size_t n);

    [[nodiscard]] const std::vector<double>& values() const noexcept {
        return values_;
    }
    [[nodiscard]] std::size_t size() const noexcept {
        return values_.size();
    }
    [[nodiscard]] double operator[](std::size_t i) const noexcept {
        return values_[i];
    }
    [[nodiscard]] double front() const noexcept {
        return values_.front();
    }
    [[nodiscard]] double back() const noexcept {
        return values_.back();
    }
    [[nodiscard]] double step() const noexcept {
        return step_;
    }

private:
    std::vector<double> values_;
    double step_ = 0.0;
};

// Domain of the E_n prior and grid (D32). The default is the toy's [10, 20]
// MeV; the source of the Geant4 runs stays inside it (E_n <= 19.5 MeV).
struct EnergyDomain {
    double En_min = 10.0 * units::MeV;
    double En_max = 20.0 * units::MeV;
    std::size_t points = kEnergyGridPoints;
};

[[nodiscard]] Grid energyGrid(const EnergyDomain& domain);

// The same lower edge and grid step as `base`, with the upper edge moved to
// En_max: the prior-robustness check widens the domain without changing the
// resolution of the grid. The point count is rounded to the nearest integer,
// so the step matches only to within the rounding of (En_max - En_min) / step.
[[nodiscard]] EnergyDomain sameStepDomain(const EnergyDomain& base, double En_max);

// theta_p on [0, pi/2]: a recoil proton never goes backwards.
[[nodiscard]] Grid recoilAngleGrid(std::size_t points = kRecoilAngleGridPoints);

} // namespace riptide::stats
