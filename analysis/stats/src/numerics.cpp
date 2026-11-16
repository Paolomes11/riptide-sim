#include "riptide/stats/numerics.hpp"

#include <cmath>
#include <limits>

namespace riptide::stats {

bool ieeeSemanticsPreserved() noexcept {
    // volatile stops the compiler from folding the checks at compile time, so
    // they observe the code generated for this translation unit.
    volatile double nan = std::numeric_limits<double>::quiet_NaN();
    volatile double inf = std::numeric_limits<double>::infinity();
    const double n = nan;
    const double i = inf;
    return std::isnan(n) && std::isinf(i) && !std::isfinite(i) && (-i < 0.0) && std::isinf(i - 1.0);
}

} // namespace riptide::stats
