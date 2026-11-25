#include "riptide/stats/direction_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

#include "riptide/common/units.hpp"

namespace riptide::stats {

namespace {

using geometry::Direction;

constexpr double kGoldenAngleFactor = 3.0;
constexpr double kHalf = 0.5;
constexpr double kStepsOfMargin = 2.0;

double goldenAngle() noexcept {
    return units::pi * (kGoldenAngleFactor - std::sqrt(5.0));
}

// Angle in [0, pi] from the clipped dot product, as the toy: the cap radius
// only needs the coarse angles, where acos is accurate enough.
double coarseAngle(const Direction& a, const Direction& b) noexcept {
    return std::acos(std::clamp(geometry::dot(a, b), -1.0, 1.0));
}

} // namespace

std::vector<Direction> sphereGrid(std::size_t n) {
    if (n == 0) {
        throw std::invalid_argument("sphereGrid: no pixels");
    }
    const double golden = goldenAngle();
    const auto count = static_cast<double>(n);
    std::vector<Direction> grid(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto index = static_cast<double>(i);
        const double cos_theta = 1.0 - 2.0 * (index + kHalf) / count;
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        const double phi = std::fmod(index * golden, 2.0 * units::pi);
        grid[i] = {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
    }
    return grid;
}

std::vector<Direction> directionCapGrid(const Direction& centre, double radius, std::size_t n) {
    if (n == 0) {
        throw std::invalid_argument("directionCapGrid: no pixels");
    }
    if (!(radius > 0.0 && radius <= units::pi)) {
        throw std::invalid_argument("directionCapGrid: radius outside (0, pi]");
    }
    const double golden = goldenAngle();
    const double one_minus_cos_radius = 1.0 - std::cos(radius);
    const auto count = static_cast<double>(n);
    std::vector<Direction> grid(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto index = static_cast<double>(i);
        const double cos_alpha = 1.0 - one_minus_cos_radius * (index + kHalf) / count;
        grid[i] = geometry::directionAroundAxis(centre, std::acos(cos_alpha), index * golden);
    }
    return grid;
}

double directionCapRadius(const std::vector<double>& log_posterior,
                          const std::vector<Direction>& grid, double delta_log) {
    if (log_posterior.size() != grid.size()) {
        throw std::invalid_argument("directionCapRadius: posterior and grid sizes differ");
    }
    if (grid.size() < 2) {
        throw std::invalid_argument("directionCapRadius: fewer than two pixels");
    }
    if (!(delta_log > 0.0) || !std::isfinite(delta_log)) {
        throw std::invalid_argument("directionCapRadius: delta_log must be positive and finite");
    }
    if (std::any_of(log_posterior.begin(), log_posterior.end(),
                    [](double v) { return std::isnan(v); })) {
        throw std::invalid_argument("directionCapRadius: NaN in the posterior");
    }
    const auto best = std::max_element(log_posterior.begin(), log_posterior.end());
    if (!std::isfinite(*best)) {
        throw std::invalid_argument("directionCapRadius: maximum is not finite");
    }
    const auto map_index = static_cast<std::size_t>(std::distance(log_posterior.begin(), best));
    const double threshold = *best - delta_log;
    const Direction& map = grid[map_index];
    double step = std::numeric_limits<double>::infinity();
    double kept = 0.0;
    for (std::size_t i = 0; i < grid.size(); ++i) {
        const double angle = coarseAngle(grid[i], map);
        if (i != map_index) {
            step = std::min(step, angle);
        }
        if (log_posterior[i] > threshold) {
            kept = std::max(kept, angle);
        }
    }
    return std::min(kept + kStepsOfMargin * step, units::pi);
}

} // namespace riptide::stats
