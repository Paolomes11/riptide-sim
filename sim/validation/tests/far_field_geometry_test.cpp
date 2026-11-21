#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>

#include "riptide/sim/beam_geometry.hpp"
#include "riptide/sim/run_config.hpp"
#include "riptide/sim/source_sampling.hpp"

// The far-field source traced through the block with the chord code of the
// checks: the fraction of rays that cross the block and the mean chord of
// those rays have closed forms, so the sampling and the geometry check each
// other.

namespace sim = riptide::sim;
namespace checks = riptide::sim::checks;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::sim::Vec3;

namespace {

constexpr std::uint64_t kSeed = 20260907; // seed of the project
constexpr double kSigmas = 5.0;
constexpr int kSamples = 200000;

const Vec3 kHalfBlock{30.0, 10.0, 10.0};            // the 60 x 20 x 20 mm block of D27
constexpr double kBlockVolume = 60.0 * 20.0 * 20.0; // [mm^3]
constexpr double kBlockSurface = 2.0 * (60.0 * 20.0 + 60.0 * 20.0 + 20.0 * 20.0); // [mm^2]

struct ChordStats {
    int hits = 0;
    double sum = 0.0;
    double sum_sq = 0.0;
    double min = 0.0;
    double max = 0.0;
};

ChordStats traceSource(const sim::SourceConfig& source) {
    std::mt19937_64 engine{kSeed};
    std::uniform_real_distribution<double> unit{0.0, 1.0};
    auto uniform = [&]() { return unit(engine); };
    const auto disk = sim::farFieldDisk(source, kHalfBlock);
    ChordStats s;
    for (int i = 0; i < kSamples; ++i) {
        const auto ray = sim::sampleFarField(source, disk, uniform);
        const auto chord = checks::rayBoxChord(ray.position, ray.direction, kHalfBlock);
        if (!chord) {
            continue;
        }
        const double length = chord->length();
        s.min = s.hits == 0 ? length : std::min(s.min, length);
        s.max = s.hits == 0 ? length : std::max(s.max, length);
        ++s.hits;
        s.sum += length;
        s.sum_sq += length * length;
    }
    return s;
}

sim::SourceConfig farField(bool isotropic, const Vec3& direction) {
    sim::SourceConfig source;
    source.beam = sim::BeamMode::FarField;
    source.isotropic = isotropic;
    source.direction = direction;
    source.margin = 10.0;
    return source;
}

// Binomial window around the expected hit fraction.
void checkHitFraction(const ChordStats& s, double expected) {
    const double sigma = std::sqrt(expected * (1.0 - expected) / kSamples);
    CHECK_THAT(static_cast<double>(s.hits) / kSamples, WithinAbs(expected, kSigmas * sigma));
}

} // namespace

TEST_CASE("far-field along +X: every hit crosses the full 60 mm", "[far_field]") {
    // Prediction: hit fraction 400 / (pi 200) = 0.6366; the shadow is the
    // 20 x 20 face, so every hit has the chord 60 mm.
    const auto source = farField(false, {1.0, 0.0, 0.0});
    const auto s = traceSource(source);
    checkHitFraction(s, sim::farFieldAcceptance(source, kHalfBlock));
    CHECK_THAT(s.min, WithinRel(60.0, 1e-12));
    CHECK_THAT(s.max, WithinRel(60.0, 1e-12));
}

TEST_CASE("far-field oblique: hit fraction and mean chord V / A", "[far_field]") {
    // Prediction: (0.6, 0, 0.8) has projected area 1200 and shadow radius
    // sqrt(1000): hit fraction 1200 / (1000 pi) = 0.3820; for a parallel beam
    // the mean chord of the hits is the volume over the projected area, 20 mm.
    const auto source = farField(false, {0.6, 0.0, 0.8});
    const auto s = traceSource(source);
    checkHitFraction(s, 1200.0 / (1000.0 * riptide::units::pi));
    const double mean = s.sum / s.hits;
    const double sigma = std::sqrt((s.sum_sq / s.hits - mean * mean) / s.hits);
    CHECK_THAT(mean, WithinAbs(kBlockVolume / 1200.0, kSigmas * sigma));
}

TEST_CASE("far-field isotropic: hit fraction S / 4 and Cauchy mean chord 4V / S", "[far_field]") {
    // Prediction: hit fraction (S / 4) / (pi R_d^2) = 1400 / (pi 1100) =
    // 0.4051; mean chord of the hits 4 V / S = 96000 / 5600 = 17.14 mm, the
    // Cauchy formula for a uniform isotropic flux on a convex body. A bias in
    // the direction weights would move this mean.
    const auto source = farField(true, {0.0, 0.0, 1.0});
    const auto s = traceSource(source);
    checkHitFraction(s, kBlockSurface / 4.0 / (riptide::units::pi * 1100.0));
    const double mean = s.sum / s.hits;
    const double sigma = std::sqrt((s.sum_sq / s.hits - mean * mean) / s.hits);
    CHECK_THAT(mean, WithinAbs(4.0 * kBlockVolume / kBlockSurface, kSigmas * sigma));
    // No chord exceeds the full diagonal of the block.
    CHECK(s.max <= 2.0 * std::sqrt(1100.0));
}
