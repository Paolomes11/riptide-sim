#pragma once

// Random numbers whose sequence is fixed by the code, not by the standard
// library: the analysis records its seed, and the same seed must give the same
// smearing with any compiler.
//
// The C++ standard fixes the output of std::mt19937_64 for a given seed, but
// not the algorithm of std::normal_distribution or std::uniform_real_distribution,
// which differ between libstdc++ and libc++. Hence the uniform and the normal
// variates are built here from the raw engine output. The normal variate still
// goes through std::log, std::sqrt and std::cos: sqrt is exactly rounded by
// IEEE 754, log and cos are not, so the last bit may differ between maths
// libraries, never the sequence of engine draws.
//
// Engines are always passed by reference: no function here owns hidden state.

#include <cmath>
#include <cstdint>
#include <random>

#include "riptide/common/units.hpp"

namespace riptide::random {

using Engine = std::mt19937_64;

// SplitMix64 (Steele, Lea and Flood, "Fast splittable pseudorandom number
// generators", OOPSLA 2014), in the form used to seed the xoshiro generators:
// the Weyl increment is the odd integer nearest 2^64 / golden ratio, and the
// two multipliers are the published output-mixing constants.
inline constexpr std::uint64_t kSplitMix64Increment = 0x9e3779b97f4a7c15ULL;
inline constexpr std::uint64_t kSplitMix64Mul1 = 0xbf58476d1ce4e5b9ULL;
inline constexpr std::uint64_t kSplitMix64Mul2 = 0x94d049bb133111ebULL;
inline constexpr unsigned kSplitMix64Shift1 = 30;
inline constexpr unsigned kSplitMix64Shift2 = 27;
inline constexpr unsigned kSplitMix64Shift3 = 31;

// Advances `state` and returns the next output. Both the step and the output
// mix are bijections of 64-bit words, so distinct states give distinct outputs.
[[nodiscard]] constexpr std::uint64_t splitMix64(std::uint64_t& state) noexcept {
    state += kSplitMix64Increment;
    std::uint64_t mixed = state;
    mixed = (mixed ^ (mixed >> kSplitMix64Shift1)) * kSplitMix64Mul1;
    mixed = (mixed ^ (mixed >> kSplitMix64Shift2)) * kSplitMix64Mul2;
    return mixed ^ (mixed >> kSplitMix64Shift3);
}

// Seed of the engine of one event. It depends only on (seed, event_id), not on
// which events come before, so a different selection or reading order leaves
// the draws of every kept event unchanged. For a fixed seed the map is
// injective in event_id: the XOR and the SplitMix64 step are bijections.
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters): (seed, event_id) as in the docs
[[nodiscard]] constexpr std::uint64_t perEventSeed(std::uint64_t seed,
                                                   std::uint64_t event_id) noexcept {
    std::uint64_t state = seed;
    state = splitMix64(state) ^ event_id;
    return splitMix64(state);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters): as perEventSeed
[[nodiscard]] inline Engine perEventEngine(std::uint64_t seed, std::uint64_t event_id) {
    return Engine{perEventSeed(seed, event_id)};
}

// A double has a 53-bit significand: the top 53 bits of a draw, scaled by
// 2^-53, give every multiple of 2^-53 in [0, 1) with equal probability.
inline constexpr unsigned kUniformBits = 53;
inline constexpr unsigned kEngineBits = 64;
inline constexpr double kUniformScale = 0x1.0p-53;

// Uniform on [0, 1).
[[nodiscard]] inline double uniform53(Engine& engine) {
    return static_cast<double>(engine() >> (kEngineBits - kUniformBits)) * kUniformScale;
}

// Standard normal by the Box-Muller transform, cosine branch only: one
// variate per call, two engine draws, so there is no cached second variate
// that would make the result depend on the call history. 1 - u lies in (0, 1],
// so the logarithm is finite.
[[nodiscard]] inline double standardNormal(Engine& engine) {
    const double u_radius = 1.0 - uniform53(engine);
    const double u_angle = uniform53(engine);
    return std::sqrt(-2.0 * std::log(u_radius)) * std::cos(2.0 * units::pi * u_angle);
}

} // namespace riptide::random
