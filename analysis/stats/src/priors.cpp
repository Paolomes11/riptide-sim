#include "riptide/stats/priors.hpp"

#include <algorithm>
#include <cmath>

#include "riptide/common/units.hpp"

namespace riptide::stats {

std::vector<double> flatEnergyLogPrior(const Grid& En_grid) {
    const double log_density = -std::log(En_grid.back() - En_grid.front());
    return std::vector<double>(En_grid.size(), log_density);
}

std::vector<double> physicalRecoilAngleLogPrior(const Grid& theta_grid) {
    constexpr double kThetaMax = 0.5 * units::pi;
    const double half = 0.5 * theta_grid.step();
    std::vector<double> log_mass(theta_grid.size());
    for (std::size_t j = 0; j < theta_grid.size(); ++j) {
        const double lo = std::clamp(theta_grid[j] - half, 0.0, kThetaMax);
        const double hi = std::clamp(theta_grid[j] + half, 0.0, kThetaMax);
        log_mass[j] = std::log(0.5 * (std::cos(2.0 * lo) - std::cos(2.0 * hi)));
    }
    return log_mass;
}

std::vector<double> flatRecoilAngleLogPrior(const Grid& theta_grid) {
    const double log_mass = -std::log(static_cast<double>(theta_grid.size()));
    return std::vector<double>(theta_grid.size(), log_mass);
}

} // namespace riptide::stats
