#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <set>
#include <vector>

#include "riptide/common/random.hpp"

namespace {

using riptide::random::Engine;
using riptide::random::perEventEngine;
using riptide::random::perEventSeed;
using riptide::random::splitMix64;
using riptide::random::standardNormal;
using riptide::random::uniform53;

constexpr std::uint64_t kProjectSeed = 20260907;
constexpr std::size_t kSamples = 1000000;
// Bounds at 4 standard errors: a correct generator fails one of these checks
// with probability ~6e-5 each, and the sequence is fixed, so a pass is stable.
constexpr double kSigmas = 4.0;
// Moments of the standard normal and the variances of their sample estimators:
// Var(mean) = 1/N, Var(s^2) = 2/N, Var(excess kurtosis) = 24/N.
constexpr double kVarianceEstimatorFactor = 2.0;
constexpr double kKurtosisEstimatorFactor = 24.0;
constexpr double kNormalKurtosis = 3.0;

std::vector<double> normalSample(std::uint64_t seed, std::size_t n) {
    Engine engine{seed};
    std::vector<double> sample(n);
    for (auto& x : sample) {
        x = standardNormal(engine);
    }
    return sample;
}

} // namespace

TEST_CASE("splitMix64 reproduces the published reference sequences", "[random]") {
    // Outputs of the reference splitmix64.c for seeds 0 and 1234567, as also
    // used by the test vectors of the xoshiro seeders.
    std::uint64_t state = 0;
    const std::array<std::uint64_t, 3> from_zero{0xe220a8397b1dcdafULL, 0x6e789e6aa1b965f4ULL,
                                                 0x06c45d188009454fULL};
    for (const auto expected : from_zero) {
        CHECK(splitMix64(state) == expected);
    }
    state = 1234567;
    const std::array<std::uint64_t, 5> from_1234567{6457827717110365317ULL, 3203168211198807973ULL,
                                                    9817491932198370423ULL, 4593380528125082431ULL,
                                                    16408922859458223821ULL};
    for (const auto expected : from_1234567) {
        CHECK(splitMix64(state) == expected);
    }
}

TEST_CASE("per-event seeds are distinct and depend only on the seed and event id", "[random]") {
    constexpr std::uint64_t kEvents = 100000;
    std::set<std::uint64_t> seeds;
    for (std::uint64_t id = 0; id < kEvents; ++id) {
        seeds.insert(perEventSeed(kProjectSeed, id));
    }
    CHECK(seeds.size() == kEvents);
    CHECK(perEventSeed(kProjectSeed, 42) == perEventSeed(kProjectSeed, 42));
    CHECK(perEventSeed(kProjectSeed, 42) != perEventSeed(kProjectSeed + 1, 42));

    Engine first = perEventEngine(kProjectSeed, 7);
    Engine again = perEventEngine(kProjectSeed, 7);
    for (int i = 0; i < 10; ++i) {
        CHECK(first() == again());
    }
}

TEST_CASE("uniform53 stays in the unit interval and has mean one half", "[random]") {
    Engine engine{kProjectSeed};
    double sum = 0.0;
    bool in_range = true;
    for (std::size_t i = 0; i < kSamples; ++i) {
        const double u = uniform53(engine);
        in_range = in_range && u >= 0.0 && u < 1.0;
        sum += u;
    }
    CHECK(in_range);
    // Var(U) = 1/12.
    const double mean_tolerance = kSigmas * std::sqrt(1.0 / 12.0 / static_cast<double>(kSamples));
    CHECK_THAT(sum / static_cast<double>(kSamples),
               Catch::Matchers::WithinAbs(0.5, mean_tolerance));
}

TEST_CASE("standardNormal has the moments of the standard normal", "[random]") {
    const auto sample = normalSample(kProjectSeed, kSamples);
    const auto n = static_cast<double>(kSamples);
    double sum = 0.0;
    for (const double x : sample) {
        sum += x;
    }
    const double mean = sum / n;
    double m2 = 0.0;
    double m4 = 0.0;
    for (const double x : sample) {
        const double d2 = (x - mean) * (x - mean);
        m2 += d2;
        m4 += d2 * d2;
    }
    m2 /= n;
    m4 /= n;

    CHECK_THAT(mean, Catch::Matchers::WithinAbs(0.0, kSigmas / std::sqrt(n)));
    CHECK_THAT(m2,
               Catch::Matchers::WithinAbs(1.0, kSigmas * std::sqrt(kVarianceEstimatorFactor / n)));
    CHECK_THAT(m4 / (m2 * m2),
               Catch::Matchers::WithinAbs(kNormalKurtosis,
                                          kSigmas * std::sqrt(kKurtosisEstimatorFactor / n)));
}

TEST_CASE("standardNormal is reproducible and seed dependent", "[random]") {
    constexpr std::size_t kShort = 1000;
    const auto a = normalSample(kProjectSeed, kShort);
    const auto b = normalSample(kProjectSeed, kShort);
    const auto c = normalSample(kProjectSeed + 1, kShort);
    bool identical = true;
    std::size_t differ = 0;
    for (std::size_t i = 0; i < kShort; ++i) {
        identical = identical && std::memcmp(&a[i], &b[i], sizeof(double)) == 0;
        differ += std::memcmp(&a[i], &c[i], sizeof(double)) != 0 ? 1U : 0U;
    }
    CHECK(identical);
    CHECK(differ == kShort);
}
