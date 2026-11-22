#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdint>
#include <random>

#include "riptide/sim/run_config.hpp"
#include "riptide/sim/source_sampling.hpp"

namespace sim = riptide::sim;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::sim::Vec3;

namespace {

constexpr std::uint64_t kSeed = 20260907; // seed of the project
constexpr double kTolerance = 1e-12;
// Width of the statistical acceptance windows, in standard deviations.
constexpr double kSigmas = 5.0;
constexpr int kSamples = 200000;

const Vec3 kHalfBlock{30.0, 10.0, 10.0}; // the 60 x 20 x 20 mm block of D27

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(const Vec3& a) {
    return std::sqrt(dot(a, a));
}

sim::SourceConfig farField(bool isotropic, const Vec3& direction) {
    sim::SourceConfig source;
    source.beam = sim::BeamMode::FarField;
    source.isotropic = isotropic;
    source.direction = direction;
    source.margin = 10.0;
    return source;
}

// Uniform doubles in [0, 1) from a seeded standard engine, the same contract
// as G4UniformRand in the generator.
struct SeededUniform {
    std::mt19937_64 engine{kSeed};
    std::uniform_real_distribution<double> unit{0.0, 1.0};
    double operator()() {
        return unit(engine);
    }
};

} // namespace

TEST_CASE("orthonormal pair completes any unit vector", "[source_sampling]") {
    const double s = 1.0 / std::sqrt(3.0);
    // -z and directions next to it are where the naive constructions fail.
    for (const Vec3& n : {Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, -1.0}, Vec3{1.0, 0.0, 0.0},
                          Vec3{s, -s, s}, Vec3{1e-9, 0.0, -std::sqrt(1.0 - 1e-18)}}) {
        const auto [e1, e2] = sim::orthonormalPair(n);
        CHECK_THAT(norm(e1), WithinAbs(1.0, kTolerance));
        CHECK_THAT(norm(e2), WithinAbs(1.0, kTolerance));
        CHECK_THAT(dot(e1, e2), WithinAbs(0.0, kTolerance));
        CHECK_THAT(dot(e1, n), WithinAbs(0.0, kTolerance));
        CHECK_THAT(dot(e2, n), WithinAbs(0.0, kTolerance));
        // Right-handed: e1 x e2 = n.
        CHECK_THAT(e1.y * e2.z - e1.z * e2.y, WithinAbs(n.x, kTolerance));
        CHECK_THAT(e1.z * e2.x - e1.x * e2.z, WithinAbs(n.y, kTolerance));
        CHECK_THAT(e1.x * e2.y - e1.y * e2.x, WithinAbs(n.z, kTolerance));
    }
}

TEST_CASE("isotropic directions are unit vectors with the moments of the sphere",
          "[source_sampling]") {
    // Prediction: on the unit sphere <Omega_i> = 0 with variance 1/3, and
    // <Omega_i^2> = 1/3 with variance <Omega_i^4> - 1/9 = 1/5 - 1/9 = 4/45.
    SeededUniform uniform;
    double sum[3] = {0.0, 0.0, 0.0};
    double sum_sq[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < kSamples; ++i) {
        const double u1 = uniform();
        const double u2 = uniform();
        const Vec3 d = sim::isotropicDirection(u1, u2);
        REQUIRE_THAT(norm(d), WithinAbs(1.0, kTolerance));
        const double c[3] = {d.x, d.y, d.z};
        for (int k = 0; k < 3; ++k) {
            sum[k] += c[k];
            sum_sq[k] += c[k] * c[k];
        }
    }
    const double n = kSamples;
    const double sigma_mean = std::sqrt(1.0 / 3.0 / n);
    const double sigma_sq = std::sqrt(4.0 / 45.0 / n);
    for (int k = 0; k < 3; ++k) {
        CHECK_THAT(sum[k] / n, WithinAbs(0.0, kSigmas * sigma_mean));
        CHECK_THAT(sum_sq[k] / n, WithinAbs(1.0 / 3.0, kSigmas * sigma_sq));
    }
    // The extremes of the uniforms map to the poles.
    CHECK_THAT(sim::isotropicDirection(0.0, 0.0).z, WithinAbs(1.0, kTolerance));
}

TEST_CASE("projected area and shadow radius of the block", "[source_sampling]") {
    // Along +X the shadow is the 20 x 20 face, held by a disk through its
    // corners: radius sqrt(10^2 + 10^2).
    const Vec3 x{1.0, 0.0, 0.0};
    CHECK_THAT(sim::projectedBoxArea(x, kHalfBlock), WithinRel(400.0, kTolerance));
    CHECK_THAT(sim::shadowRadius(x, kHalfBlock), WithinRel(std::sqrt(200.0), kTolerance));
    CHECK_THAT(sim::projectedBoxArea({0.0, 1.0, 0.0}, kHalfBlock), WithinRel(1200.0, kTolerance));
    // Oblique (0.6, 0, 0.8): area 4 (0.6 * 100 + 0.8 * 300) = 1200; the corner
    // projections are 18 +- 8, the smallest is 10, radius^2 = 1100 - 100.
    const Vec3 oblique{0.6, 0.0, 0.8};
    CHECK_THAT(sim::projectedBoxArea(oblique, kHalfBlock), WithinRel(1200.0, kTolerance));
    CHECK_THAT(sim::shadowRadius(oblique, kHalfBlock), WithinRel(std::sqrt(1000.0), kTolerance));
    CHECK_THAT(sim::halfDiagonal(kHalfBlock), WithinRel(std::sqrt(1100.0), kTolerance));

    // For any direction every corner lies within the shadow radius of the line.
    SeededUniform uniform;
    for (int i = 0; i < 1000; ++i) {
        const double u1 = uniform();
        const double u2 = uniform();
        const Vec3 d = sim::isotropicDirection(u1, u2);
        const double radius = sim::shadowRadius(d, kHalfBlock);
        for (const double sx : {-1.0, 1.0}) {
            for (const double sy : {-1.0, 1.0}) {
                for (const double sz : {-1.0, 1.0}) {
                    const Vec3 c{sx * kHalfBlock.x, sy * kHalfBlock.y, sz * kHalfBlock.z};
                    const double along = dot(c, d);
                    const double off_line_sq = dot(c, c) - along * along;
                    CHECK(off_line_sq <= radius * radius * (1.0 + kTolerance));
                }
            }
        }
    }
}

TEST_CASE("far-field disk and acceptance", "[source_sampling]") {
    const double diagonal = std::sqrt(1100.0);
    const auto iso = farField(true, {0.0, 0.0, 1.0});
    const auto disk_iso = sim::farFieldDisk(iso, kHalfBlock);
    CHECK_THAT(disk_iso.radius, WithinRel(diagonal, kTolerance));
    CHECK_THAT(disk_iso.distance, WithinRel(diagonal + 10.0, kTolerance));
    // Mean shadow S / 4 = 1400 over the disk area pi 1100.
    CHECK_THAT(sim::farFieldAcceptance(iso, kHalfBlock),
               WithinRel(1400.0 / (riptide::units::pi * 1100.0), kTolerance));

    const auto fixed = farField(false, {1.0, 0.0, 0.0});
    CHECK_THAT(sim::farFieldDisk(fixed, kHalfBlock).radius,
               WithinRel(std::sqrt(200.0), kTolerance));
    CHECK_THAT(sim::farFieldAcceptance(fixed, kHalfBlock),
               WithinRel(400.0 / (riptide::units::pi * 200.0), kTolerance));

    CHECK_THROWS_AS(sim::farFieldDisk(sim::SourceConfig{}, kHalfBlock), sim::ConfigError);
}

TEST_CASE("far-field start points: on the disk, uniform in area, outside the block",
          "[source_sampling]") {
    // Prediction: the start point minus the disk centre is perpendicular to
    // the direction and within the radius; (rho / r)^2 is uniform, mean 1/2
    // and variance 1/12. The disk lies beyond the circumscribed sphere, so no
    // start point is inside the block.
    const auto source = farField(true, {0.0, 0.0, 1.0});
    const auto disk = sim::farFieldDisk(source, kHalfBlock);
    SeededUniform uniform;
    double sum_area_fraction = 0.0;
    for (int i = 0; i < kSamples; ++i) {
        const auto ray = sim::sampleFarField(source, disk, uniform);
        const Vec3 offset{ray.position.x + disk.distance * ray.direction.x,
                          ray.position.y + disk.distance * ray.direction.y,
                          ray.position.z + disk.distance * ray.direction.z};
        REQUIRE_THAT(dot(offset, ray.direction), WithinAbs(0.0, 1e-10));
        const double rho_sq = dot(offset, offset);
        REQUIRE(rho_sq <= disk.radius * disk.radius * (1.0 + kTolerance));
        REQUIRE(norm(ray.position) > disk.distance * (1.0 - kTolerance));
        sum_area_fraction += rho_sq / (disk.radius * disk.radius);
    }
    CHECK_THAT(sum_area_fraction / kSamples,
               WithinAbs(0.5, kSigmas * std::sqrt(1.0 / 12.0 / kSamples)));
}

TEST_CASE("far-field sampling draws two uniforms per random quantity", "[source_sampling]") {
    int calls = 0;
    auto counting = [&calls]() {
        ++calls;
        return 0.25;
    };
    const auto iso = farField(true, {0.0, 0.0, 1.0});
    (void)sim::sampleFarField(iso, sim::farFieldDisk(iso, kHalfBlock), counting);
    CHECK(calls == 4);

    calls = 0;
    const auto fixed = farField(false, {1.0, 0.0, 0.0});
    const auto ray = sim::sampleFarField(fixed, sim::farFieldDisk(fixed, kHalfBlock), counting);
    CHECK(calls == 2);
    CHECK_THAT(ray.direction.x, WithinAbs(1.0, 0.0));
}

TEST_CASE("source extent bounds every start point", "[source_sampling]") {
    sim::SourceConfig pencil;
    pencil.position = {3.0, -120.0, 7.0};
    CHECK_THAT(sim::sourceExtent(pencil, kHalfBlock), WithinAbs(120.0, 0.0));

    const auto source = farField(true, {0.0, 0.0, 1.0});
    const auto disk = sim::farFieldDisk(source, kHalfBlock);
    const double extent = sim::sourceExtent(source, kHalfBlock);
    SeededUniform uniform;
    for (int i = 0; i < 10000; ++i) {
        const auto ray = sim::sampleFarField(source, disk, uniform);
        REQUIRE(norm(ray.position) <= extent * (1.0 + kTolerance));
    }
}
