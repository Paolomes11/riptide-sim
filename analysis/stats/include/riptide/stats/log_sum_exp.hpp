#pragma once

// Numerically stable log(sum exp(x_i)). Marginalising theta_p sums terms whose
// logarithms reach -1e4 for events far from a grid cell; exponentiating them
// directly underflows to 0 and the log of the sum becomes -inf for every E_n.
// Shifting by the maximum keeps the largest term at exp(0) = 1.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// -inf terms are allowed (impossible configurations) and an all -inf or empty
// input gives -inf. NaN and +inf are not log-probabilities: std::invalid_argument.
[[nodiscard]] inline double logSumExp(const double* values, std::size_t n) {
    constexpr double kMinusInf = -std::numeric_limits<double>::infinity();
    double max_value = kMinusInf;
    for (std::size_t i = 0; i < n; ++i) {
        if (std::isnan(values[i]) || (std::isinf(values[i]) && values[i] > 0.0)) {
            throw std::invalid_argument("logSumExp: NaN or +inf term");
        }
        max_value = std::max(max_value, values[i]);
    }
    if (std::isinf(max_value)) {
        return kMinusInf;
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += std::exp(values[i] - max_value);
    }
    return max_value + std::log(sum);
}

// The same sum without the NaN / +inf check, for hot loops whose terms are
// finite or -inf by construction. A NaN or +inf input gives an unspecified
// result instead of an exception.
[[nodiscard]] inline double logSumExpUnchecked(const double* values, std::size_t n) noexcept {
    constexpr double kMinusInf = -std::numeric_limits<double>::infinity();
    double max_value = kMinusInf;
    for (std::size_t i = 0; i < n; ++i) {
        max_value = std::max(max_value, values[i]);
    }
    if (std::isinf(max_value)) {
        return kMinusInf;
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += std::exp(values[i] - max_value);
    }
    return max_value + std::log(sum);
}

[[nodiscard]] inline double logSumExp(const std::vector<double>& values) {
    return logSumExp(values.data(), values.size());
}

} // namespace riptide::stats
