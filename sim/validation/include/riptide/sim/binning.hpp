#pragma once

// Equal-width binning shared by the phase 1b and phase 2 checks, so both
// report histograms with the same edge convention.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace riptide::sim::checks {

// `bins` + 1 edges spanning [low, high]; the last one is exactly `high`,
// whatever the rounding of the step.
[[nodiscard]] inline std::vector<double> uniformEdges(double low, double high, std::size_t bins) {
    std::vector<double> edges(bins + 1);
    for (std::size_t i = 0; i <= bins; ++i) {
        edges.at(i) = low + (high - low) * static_cast<double>(i) / static_cast<double>(bins);
    }
    edges.back() = high;
    return edges;
}

// Bin of x among `bins` equal bins of [0, high], the upper edge in the last bin.
// The caller keeps x in [0, high]: values above it would land in the last bin.
[[nodiscard]] inline std::size_t binIndex(double x, double high, std::size_t bins) {
    const auto i = static_cast<std::size_t>(std::floor(x / high * static_cast<double>(bins)));
    return std::min(i, bins - 1);
}

} // namespace riptide::sim::checks
