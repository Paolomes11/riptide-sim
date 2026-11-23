#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <set>
#include <stdexcept>
#include <vector>

#include "riptide/common/random.hpp"

namespace {

using riptide::random::Engine;
using riptide::random::perEventEngine;
using riptide::random::perEventSeed;
using riptide::random::rejectionThreshold;
using riptide::random::splitMix64;
using riptide::random::standardNormal;
using riptide::random::uniform53;
using riptide::random::uniformBelow;

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

TEST_CASE("rejectionThreshold is 2^64 mod bound", "[random]") {
    // 2^64 = (2^2)^32 is 1 mod 3; a power of two divides 2^64; and
    // 2^64 = 2 (2^63 + 1) - 2, so the residue is 2^63 - 1, the largest
    // rejected range of any bound.
    constexpr std::uint64_t kThree = 3;
    constexpr std::uint64_t kPowerOfTwo = std::uint64_t{1} << 40U;
    constexpr std::uint64_t kHalfPlusOne = (std::uint64_t{1} << 63U) + 1;
    STATIC_CHECK(rejectionThreshold(1) == 0);
    STATIC_CHECK(rejectionThreshold(kThree) == 1);
    STATIC_CHECK(rejectionThreshold(kPowerOfTwo) == 0);
    STATIC_CHECK(rejectionThreshold(kHalfPlusOne) == (std::uint64_t{1} << 63U) - 1);
}

TEST_CASE("uniformBelow stays in range and rejects a zero bound", "[random]") {
    Engine engine{kProjectSeed};
    CHECK_THROWS_AS(uniformBelow(engine, 0), std::invalid_argument);
    constexpr std::size_t kDraws = 1000;
    constexpr std::uint64_t kHalfPlusOne = (std::uint64_t{1} << 63U) + 1;
    bool ones_give_zero = true;
    bool below_bound = true;
    for (std::size_t i = 0; i < kDraws; ++i) {
        ones_give_zero = ones_give_zero && uniformBelow(engine, 1) == 0;
        below_bound = below_bound && uniformBelow(engine, kHalfPlusOne) < kHalfPlusOne;
    }
    CHECK(ones_give_zero);
    CHECK(below_bound);
}

TEST_CASE("uniformBelow gives every value with equal frequency", "[random]") {
    // 7 does not divide 2^64, so the rejection step is exercised.
    constexpr std::uint64_t kBound = 7;
    Engine engine{kProjectSeed};
    std::array<std::size_t, kBound> counts{};
    for (std::size_t i = 0; i < kSamples; ++i) {
        ++counts.at(uniformBelow(engine, kBound));
    }
    const double n = static_cast<double>(kSamples);
    const double p = 1.0 / static_cast<double>(kBound);
    const double sigma = std::sqrt(p * (1.0 - p) / n);
    for (const std::size_t count : counts) {
        CHECK_THAT(static_cast<double>(count) / n, Catch::Matchers::WithinAbs(p, kSigmas * sigma));
    }
}

TEST_CASE("uniformBelow is reproducible and seed dependent", "[random]") {
    constexpr std::size_t kShort = 1000;
    constexpr std::uint64_t kBound = 1000003;
    Engine a{kProjectSeed};
    Engine b{kProjectSeed};
    Engine c{kProjectSeed + 1};
    bool identical = true;
    std::size_t differ = 0;
    for (std::size_t i = 0; i < kShort; ++i) {
        const std::uint64_t x = uniformBelow(a, kBound);
        identical = identical && x == uniformBelow(b, kBound);
        differ += x != uniformBelow(c, kBound) ? 1U : 0U;
    }
    CHECK(identical);
    // Two independent streams agree at a given draw with probability 1/kBound.
    CHECK(differ + 1 >= kShort);
}
