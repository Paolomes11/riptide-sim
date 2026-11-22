#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "riptide/common/event_sample.hpp"
#include "riptide/common/hp_kinematics.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"
#include "riptide/io/response.hpp"
#include "riptide/io/response_summary.hpp"
#include "riptide/io/run_info.hpp"

namespace {

using riptide::EventSample;
using riptide::io::IoError;
using riptide::io::Resolution;
using riptide::kinematics::kNeutronMass;

constexpr std::uint64_t kProjectSeed = 20260907;
constexpr double kSigmaEp = 0.10;    // MeV
constexpr double kSigmaTheta = 0.08; // rad
constexpr double kAwrH1 = 0.999167;
constexpr std::size_t kSamples = 100000;
// A correct pull is N(0, 1); 4 sigma keeps the test from flaking.
constexpr double kPullBound = 4.0;
constexpr double kOffsetMin = 5.33500e-4;
constexpr double kOffsetMax = 5.50179e-4;
constexpr double kLabAngleTolerance = 6e-6;

// HP events at 14 MeV over the forward hemisphere, with tracks in the x-z plane.
EventSample hpSample(std::size_t n) {
    EventSample sample;
    riptide::random::Engine angles{kProjectSeed};
    const double target_mass = kAwrH1 * kNeutronMass;
    constexpr double kEn = 14.0;
    // Stay below 90 degrees, where the proton energy vanishes.
    constexpr double kThetaMax = 1.5;
    for (std::size_t i = 0; i < n; ++i) {
        const double theta_p = kThetaMax * riptide::random::uniform53(angles);
        sample.event_id.push_back(static_cast<std::int32_t>(i));
        sample.truth_En.push_back(kEn);
        sample.truth_omega_n_x.push_back(0.0);
        sample.truth_omega_n_y.push_back(0.0);
        sample.truth_omega_n_z.push_back(1.0);
        sample.truth_theta_p.push_back(theta_p);
        sample.truth_Ep.push_back(riptide::kinematics::hpProtonEnergy(kEn, theta_p, target_mass));
        sample.truth_p_dir_x.push_back(std::sin(theta_p));
        sample.truth_p_dir_y.push_back(0.0);
        sample.truth_p_dir_z.push_back(std::cos(theta_p));
        sample.det_Evis.push_back(0.0);
    }
    riptide::io::convertHpEnergies(sample, target_mass);
    return sample;
}

} // namespace

TEST_CASE("the conversion summary reports the HP offset", "[response_summary]") {
    const auto sample = hpSample(1000);
    const auto doc = riptide::io::hpConversionSummary(sample);
    CHECK(doc.at("events") == 1000);
    CHECK(doc.at("events_without_offset") == 0);
    const auto& offset = doc.at("relative_offset");
    CHECK(offset.at("min").get<double>() > kOffsetMin);
    CHECK(offset.at("max").get<double>() < kOffsetMax);
    CHECK(offset.at("min").get<double>() <= offset.at("mean").get<double>());
    CHECK(offset.at("mean").get<double>() <= offset.at("max").get<double>());
    CHECK(doc.at("lab_angle_relative_max").get<double>() < kLabAngleTolerance);
}

TEST_CASE("the conversion summary needs the physical energies", "[response_summary]") {
    auto sample = hpSample(10);
    sample.truth_Ep_phys.clear();
    CHECK_THROWS_AS(riptide::io::hpConversionSummary(sample), IoError);
    CHECK_THROWS_AS(riptide::io::hpConversionSummary(EventSample{}), IoError);
}

TEST_CASE("smearing pulls are standard for the applied widths", "[response_summary]") {
    auto sample = hpSample(kSamples);
    const Resolution resolution{kSigmaEp, kSigmaTheta};
    riptide::io::applyResolution(sample, resolution, kProjectSeed);
    const auto doc = riptide::io::resolutionSummary(sample, resolution);
    CHECK(doc.at("events") == kSamples);
    for (const char* key : {"Ep_MeV", "theta_p_rad"}) {
        INFO(key);
        CHECK(std::abs(doc.at(key).at("pull_mean").get<double>()) < kPullBound);
        CHECK(std::abs(doc.at(key).at("pull_variance").get<double>()) < kPullBound);
    }
    CHECK_THAT(doc.at("Ep_MeV").at("std").get<double>(),
               Catch::Matchers::WithinRel(kSigmaEp, kPullBound / std::sqrt(2.0 * kSamples)));
    const auto& direction = doc.at("p_dir");
    CHECK(std::abs(direction.at("pull_mean").get<double>()) < kPullBound);
    CHECK(direction.at("norm_deviation_max").get<double>() < 1e-12);
}

TEST_CASE("a wrong width shows up in the pulls", "[response_summary]") {
    auto sample = hpSample(kSamples);
    riptide::io::applyResolution(sample, Resolution{kSigmaEp, kSigmaTheta}, kProjectSeed);
    // Judged against widths 10% too small, the variances are far off.
    constexpr double kWrongScale = 0.9;
    const auto doc = riptide::io::resolutionSummary(
        sample, Resolution{kWrongScale * kSigmaEp, kWrongScale * kSigmaTheta});
    CHECK(doc.at("Ep_MeV").at("pull_variance").get<double>() > 10.0 * kPullBound);
    CHECK(doc.at("theta_p_rad").at("pull_variance").get<double>() > 10.0 * kPullBound);
    CHECK(doc.at("p_dir").at("pull_mean").get<double>() > 10.0 * kPullBound);
}

TEST_CASE("zero widths give zero residuals and no pulls", "[response_summary]") {
    auto sample = hpSample(100);
    const Resolution none{0.0, 0.0};
    riptide::io::applyResolution(sample, none, kProjectSeed);
    const auto doc = riptide::io::resolutionSummary(sample, none);
    CHECK_THAT(doc.at("Ep_MeV").at("max_abs").get<double>(), Catch::Matchers::WithinAbs(0.0, 0.0));
    CHECK_THAT(doc.at("theta_p_rad").at("max_abs").get<double>(),
               Catch::Matchers::WithinAbs(0.0, 0.0));
    CHECK_THAT(doc.at("p_dir").at("angle_max").get<double>(), Catch::Matchers::WithinAbs(0.0, 0.0));
    CHECK(doc.at("Ep_MeV").at("pull_mean").is_null());
    CHECK(doc.at("p_dir").at("pull_mean").is_null());
}

TEST_CASE("the resolution summary needs the reco columns", "[response_summary]") {
    auto sample = hpSample(10);
    CHECK_THROWS_AS(riptide::io::resolutionSummary(sample, Resolution{kSigmaEp, kSigmaTheta}),
                    IoError);
    riptide::io::applyResolution(sample, Resolution{kSigmaEp, kSigmaTheta}, kProjectSeed);
    auto single = sample;
    single.forEachDoubleColumn([](std::vector<double>& column) { column.resize(1); });
    single.event_id.resize(1);
    CHECK_THROWS_AS(riptide::io::resolutionSummary(single, Resolution{kSigmaEp, kSigmaTheta}),
                    IoError);
}

TEST_CASE("events without a track leave the direction pull defined", "[response_summary]") {
    constexpr std::size_t kWithoutTrack = 50;
    auto sample = hpSample(kSamples);
    // The first events as the simulation writes them when no proton was made.
    for (std::size_t i = 0; i < kWithoutTrack; ++i) {
        sample.truth_Ep[i] = 0.0;
        sample.truth_Ep_phys[i] = 0.0;
        sample.truth_theta_p[i] = 0.0;
        sample.truth_p_dir_x[i] = 0.0;
        sample.truth_p_dir_y[i] = 0.0;
        sample.truth_p_dir_z[i] = 0.0;
    }
    const Resolution resolution{kSigmaEp, kSigmaTheta};
    riptide::io::applyResolution(sample, resolution, kProjectSeed);
    const auto doc = riptide::io::resolutionSummary(sample, resolution);
    const auto& dir = doc.at("p_dir");
    CHECK(dir.at("events_without_track") == kWithoutTrack);
    REQUIRE(dir.at("pull_mean").is_number());
    CHECK(std::abs(dir.at("pull_mean").get<double>()) < kPullBound);
    CHECK(riptide::io::hpConversionSummary(sample).at("events_without_offset") == kWithoutTrack);
}

TEST_CASE("a non-finite residual is refused, not hidden by a maximum", "[response_summary]") {
    const Resolution resolution{kSigmaEp, kSigmaTheta};
    auto sample = hpSample(10);
    riptide::io::applyResolution(sample, resolution, kProjectSeed);
    auto bad_dir = sample;
    bad_dir.reco_p_dir_x[3] = NAN;
    CHECK_THROWS_AS(riptide::io::resolutionSummary(bad_dir, resolution), IoError);
    auto bad_Ep = sample;
    bad_Ep.reco_Ep[3] = NAN;
    CHECK_THROWS_AS(riptide::io::resolutionSummary(bad_Ep, resolution), IoError);
    auto bad_truth = sample;
    bad_truth.truth_p_dir_z[3] = 0.5;
    CHECK_THROWS_AS(riptide::io::resolutionSummary(bad_truth, resolution), IoError);
}
