#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "riptide/common/event_sample.hpp"
#include "riptide/common/hp_kinematics.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"
#include "riptide/io/response.hpp"

namespace {

using riptide::EventSample;
using riptide::io::Direction;
using riptide::io::Resolution;
using riptide::kinematics::kNeutronMass;
using riptide::kinematics::kProtonMass;
using riptide::random::Engine;

constexpr std::uint64_t kProjectSeed = 20260907;
// Resolutions of the toy model.
constexpr double kSigmaEp = 0.10;    // MeV
constexpr double kSigmaTheta = 0.08; // rad
// AWR of H-1 in G4NDL 4.7.1.
constexpr double kAwrH1 = 0.999167;
constexpr double kUnitTolerance = 1e-12;
constexpr std::size_t kSamples = 200000;
constexpr double kSigmas = 4.0;

using riptide::geometry::angleBetween;
using riptide::geometry::dot;

double norm(const Direction& v) {
    return std::sqrt(dot(v, v));
}

Direction randomUnit(Engine& engine) {
    const double cos_theta = 2.0 * riptide::random::uniform53(engine) - 1.0;
    const double phi = 2.0 * riptide::units::pi * riptide::random::uniform53(engine);
    const double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
}

EventSample sampleWithIds(const std::vector<std::int32_t>& ids) {
    EventSample sample;
    Engine engine{kProjectSeed};
    for (const auto id : ids) {
        sample.event_id.push_back(id);
        const auto track = randomUnit(engine);
        sample.truth_En.push_back(14.0);
        sample.truth_Ep.push_back(5.0 + 0.001 * id);
        sample.truth_Ep_phys.push_back(5.0 + 0.001 * id);
        sample.truth_theta_p.push_back(std::acos(track.z));
        sample.truth_p_dir_x.push_back(track.x);
        sample.truth_p_dir_y.push_back(track.y);
        sample.truth_p_dir_z.push_back(track.z);
    }
    return sample;
}

bool sameBits(double a, double b) {
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

} // namespace

TEST_CASE("the perpendicular basis is a right-handed orthonormal frame", "[response]") {
    const std::array<Direction, 4> axes{Direction{1.0, 0.0, 0.0}, Direction{0.0, 1.0, 0.0},
                                        Direction{0.0, 0.0, -1.0}, Direction{0.6, 0.0, 0.8}};
    for (const auto& axis : axes) {
        const auto [e1, e2] = riptide::io::perpendicularBasis(axis);
        CHECK_THAT(norm(e1), Catch::Matchers::WithinAbs(1.0, kUnitTolerance));
        CHECK_THAT(norm(e2), Catch::Matchers::WithinAbs(1.0, kUnitTolerance));
        CHECK_THAT(dot(e1, axis), Catch::Matchers::WithinAbs(0.0, kUnitTolerance));
        CHECK_THAT(dot(e2, axis), Catch::Matchers::WithinAbs(0.0, kUnitTolerance));
        CHECK_THAT(dot(e1, e2), Catch::Matchers::WithinAbs(0.0, kUnitTolerance));
        // e1 x e2 = axis for a right-handed frame.
        const Direction e3{e1.y * e2.z - e1.z * e2.y, e1.z * e2.x - e1.x * e2.z,
                           e1.x * e2.y - e1.y * e2.x};
        CHECK_THAT(dot(e3, axis), Catch::Matchers::WithinAbs(1.0, kUnitTolerance));
    }
}

TEST_CASE("directionAroundAxis keeps the polar angle", "[response]") {
    const Direction axis{0.0, 0.6, 0.8};
    const double theta = 0.3;
    for (const double phi : {0.0, 1.0, 2.5, -2.0}) {
        const auto direction = riptide::io::directionAroundAxis(axis, theta, phi);
        CHECK_THAT(norm(direction), Catch::Matchers::WithinAbs(1.0, kUnitTolerance));
        CHECK_THAT(angleBetween(axis, direction),
                   Catch::Matchers::WithinAbs(theta, kUnitTolerance));
    }
}

TEST_CASE("smearDirection gives unit vectors at a Rayleigh angle", "[response]") {
    Engine tracks{kProjectSeed};
    Engine noise{kProjectSeed + 1};
    double sum_angle_squared = 0.0;
    std::size_t within_sigma = 0;
    bool unit = true;
    for (std::size_t i = 0; i < kSamples; ++i) {
        const auto track = randomUnit(tracks);
        const auto smeared = riptide::io::smearDirection(track, kSigmaTheta, noise);
        unit = unit && std::abs(norm(smeared) - 1.0) < kUnitTolerance;
        const double angle = angleBetween(track, smeared);
        sum_angle_squared += angle * angle;
        within_sigma += angle < kSigmaTheta ? 1U : 0U;
    }
    CHECK(unit);
    const auto n = static_cast<double>(kSamples);
    // angle^2 / (2 sigma^2) is Exp(1): mean 1, standard deviation 1.
    const double mean_scaled = sum_angle_squared / n / (2.0 * kSigmaTheta * kSigmaTheta);
    CHECK_THAT(mean_scaled, Catch::Matchers::WithinAbs(1.0, kSigmas / std::sqrt(n)));
    // Rayleigh CDF at sigma: 1 - exp(-1/2).
    const double p_within = 1.0 - std::exp(-0.5);
    CHECK_THAT(
        static_cast<double>(within_sigma) / n,
        Catch::Matchers::WithinAbs(p_within, kSigmas * std::sqrt(p_within * (1.0 - p_within) / n)));
}

TEST_CASE("zero resolution leaves the truth unchanged", "[response]") {
    auto sample = sampleWithIds({0, 1, 2, 3});
    riptide::io::applyResolution(sample, Resolution{0.0, 0.0}, kProjectSeed);
    for (std::size_t i = 0; i < sample.size(); ++i) {
        CHECK(sameBits(sample.reco_Ep[i], sample.truth_Ep_phys[i]));
        CHECK(sameBits(sample.reco_theta_p[i], sample.truth_theta_p[i]));
        CHECK(sameBits(sample.reco_p_dir_x[i], sample.truth_p_dir_x[i]));
        CHECK(sameBits(sample.reco_p_dir_y[i], sample.truth_p_dir_y[i]));
        CHECK(sameBits(sample.reco_p_dir_z[i], sample.truth_p_dir_z[i]));
    }
}

TEST_CASE("the noise of an event does not depend on the other events", "[response]") {
    auto full = sampleWithIds({3, 7, 11, 19});
    // Same truth for event 7, alone and at another position.
    auto alone = sampleWithIds({3, 7, 11, 19});
    // Keep only the second event; columns not filled by the fixture stay empty.
    alone.forEachDoubleColumn([](std::vector<double>& column) {
        if (!column.empty()) {
            column = {column.at(1)};
        }
    });
    alone.event_id = {7};
    const Resolution resolution{kSigmaEp, kSigmaTheta};
    riptide::io::applyResolution(full, resolution, kProjectSeed);
    riptide::io::applyResolution(alone, resolution, kProjectSeed);
    CHECK(sameBits(full.reco_Ep[1], alone.reco_Ep[0]));
    CHECK(sameBits(full.reco_theta_p[1], alone.reco_theta_p[0]));
    CHECK(sameBits(full.reco_p_dir_x[1], alone.reco_p_dir_x[0]));
    CHECK(sameBits(full.reco_p_dir_z[1], alone.reco_p_dir_z[0]));

    auto other_seed = sampleWithIds({3, 7, 11, 19});
    riptide::io::applyResolution(other_seed, resolution, kProjectSeed + 1);
    CHECK_FALSE(sameBits(full.reco_Ep[1], other_seed.reco_Ep[1]));
}

TEST_CASE("invalid resolutions and inputs are refused", "[response]") {
    auto sample = sampleWithIds({0, 1});
    CHECK_THROWS_AS(riptide::io::applyResolution(sample, Resolution{-0.1, kSigmaTheta}, 1),
                    std::invalid_argument);
    CHECK_THROWS_AS(riptide::io::applyResolution(sample, Resolution{kSigmaEp, NAN}, 1),
                    std::invalid_argument);
    sample.event_id[1] = -1;
    CHECK_THROWS_AS(riptide::io::applyResolution(sample, Resolution{kSigmaEp, kSigmaTheta}, 1),
                    std::invalid_argument);
    sample = sampleWithIds({0, 1});
    sample.truth_theta_p.pop_back();
    CHECK_THROWS_AS(riptide::io::applyResolution(sample, Resolution{kSigmaEp, kSigmaTheta}, 1),
                    std::invalid_argument);
    CHECK_THROWS_AS(riptide::io::convertHpEnergies(sample, 0.0), std::invalid_argument);
}

TEST_CASE("with the proton mass the HP conversion is the identity", "[response]") {
    auto sample = sampleWithIds({0, 1, 2, 3, 4});
    riptide::io::convertHpEnergies(sample, kProtonMass);
    for (std::size_t i = 0; i < sample.size(); ++i) {
        CHECK_THAT(sample.truth_Ep_phys[i],
                   Catch::Matchers::WithinRel(sample.truth_Ep[i], kUnitTolerance));
    }
}

TEST_CASE("the HP conversion removes the expected offset", "[response]") {
    // Expected relative offset (Ep_hp - Ep_phys) / Ep_phys over 10-19.5 MeV,
    // from the phase-3 predictions; and the physical energy agrees with the
    // bare-proton kinematics at the same lab angle to a few 1e-6 at most.
    constexpr double kOffsetMin = 5.33500e-4;
    constexpr double kOffsetMax = 5.50179e-4;
    constexpr double kLabAngleTolerance = 6e-6;
    const double target_mass = kAwrH1 * kNeutronMass;
    EventSample sample;
    for (const double En : {10.0, 14.0, 19.5}) {
        for (const double theta_p : {0.05, 0.4, 0.8, 1.2}) {
            sample.event_id.push_back(static_cast<std::int32_t>(sample.event_id.size()));
            sample.truth_En.push_back(En);
            sample.truth_theta_p.push_back(theta_p);
            sample.truth_Ep.push_back(
                riptide::kinematics::hpProtonEnergy(En, theta_p, target_mass));
        }
    }
    riptide::io::convertHpEnergies(sample, target_mass);
    for (std::size_t i = 0; i < sample.size(); ++i) {
        const double offset =
            (sample.truth_Ep[i] - sample.truth_Ep_phys[i]) / sample.truth_Ep_phys[i];
        CHECK(offset > kOffsetMin);
        CHECK(offset < kOffsetMax);
        CHECK_THAT(sample.truth_Ep_phys[i],
                   Catch::Matchers::WithinRel(riptide::kinematics::protonEnergy(
                                                  sample.truth_En[i], sample.truth_theta_p[i]),
                                              kLabAngleTolerance));
    }
}

TEST_CASE("a zero truth direction means no track", "[response]") {
    CHECK(riptide::io::hasTrack(Direction{0.0, 0.0, 1.0}));
    CHECK(riptide::io::hasTrack(Direction{0.6, 0.0, 0.8}));
    CHECK_FALSE(riptide::io::hasTrack(Direction{0.0, 0.0, 0.0}));
    CHECK_THROWS_AS(riptide::io::hasTrack(Direction{0.5, 0.0, 0.0}), std::invalid_argument);
    CHECK_THROWS_AS(riptide::io::hasTrack(Direction{NAN, 0.0, 1.0}), std::invalid_argument);
}

TEST_CASE("an event without a track gets a zero reco direction", "[response]") {
    const Resolution resolution{kSigmaEp, kSigmaTheta};
    auto with_tracks = sampleWithIds({0, 1, 2});
    auto without = sampleWithIds({0, 1, 2});
    // Event 1 as the simulation writes it when the first n-p made no proton.
    without.truth_Ep[1] = 0.0;
    without.truth_Ep_phys[1] = 0.0;
    without.truth_theta_p[1] = 0.0;
    without.truth_p_dir_x[1] = 0.0;
    without.truth_p_dir_y[1] = 0.0;
    without.truth_p_dir_z[1] = 0.0;
    riptide::io::applyResolution(with_tracks, resolution, kProjectSeed);
    riptide::io::applyResolution(without, resolution, kProjectSeed);

    CHECK_THAT(without.reco_p_dir_x[1], Catch::Matchers::WithinAbs(0.0, 0.0));
    CHECK_THAT(without.reco_p_dir_y[1], Catch::Matchers::WithinAbs(0.0, 0.0));
    CHECK_THAT(without.reco_p_dir_z[1], Catch::Matchers::WithinAbs(0.0, 0.0));
    // Ep and theta_p are drawn before the direction: their noise is the same,
    // up to the rounding of (truth + noise) - truth.
    CHECK_THAT(without.reco_Ep[1] - without.truth_Ep_phys[1],
               Catch::Matchers::WithinAbs(with_tracks.reco_Ep[1] - with_tracks.truth_Ep_phys[1],
                                          kUnitTolerance));
    CHECK_THAT(without.reco_theta_p[1] - without.truth_theta_p[1],
               Catch::Matchers::WithinAbs(
                   with_tracks.reco_theta_p[1] - with_tracks.truth_theta_p[1], kUnitTolerance));
    // The other events do not notice.
    for (const std::size_t i : {std::size_t{0}, std::size_t{2}}) {
        CHECK(sameBits(without.reco_p_dir_x[i], with_tracks.reco_p_dir_x[i]));
        CHECK(sameBits(without.reco_p_dir_z[i], with_tracks.reco_p_dir_z[i]));
    }

    auto damaged = sampleWithIds({0, 1});
    damaged.truth_p_dir_x[1] = NAN;
    CHECK_THROWS_AS(riptide::io::applyResolution(damaged, resolution, kProjectSeed),
                    std::invalid_argument);
}
