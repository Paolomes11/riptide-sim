#pragma once

// Floating-point guarantees the statistics code relies on. Posterior grids and
// log-likelihood sums use -inf for impossible configurations and NaN checks to
// catch invalid input; -ffast-math (or -ffinite-math-only) lets the compiler
// assume neither exists and silently breaks both. Every translation unit of
// this library refuses to build under those flags.

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
#error                                                                                             \
    "riptide::stats requires IEEE semantics: do not build it with -ffast-math or -ffinite-math-only"
#endif

namespace riptide::stats {

// True when NaN and infinities behave as IEEE 754 prescribes at run time.
[[nodiscard]] bool ieeeSemanticsPreserved() noexcept;

} // namespace riptide::stats
