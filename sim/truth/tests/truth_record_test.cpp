#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <optional>

#include "riptide/common/units.hpp"
#include "riptide/sim/truth_record.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::schema::FirstInteraction;
using riptide::sim::classifyFirstInteraction;
using riptide::sim::HadronicKind;
using riptide::sim::HadronicStep;
using riptide::sim::kCarbonZ;
using riptide::sim::kHydrogenZ;
using riptide::sim::RecoilProton;
using riptide::sim::thetaP;
using riptide::sim::TruthRecorder;
using riptide::sim::TruthTally;
using riptide::sim::Vec3;
namespace u = riptide::units;

constexpr double kAngleTol = 1e-12;
constexpr int kOxygenZ = 8;
constexpr double kEn = 14.0 * u::MeV;
const Vec3 kBeam{0.0, 0.0, 1.0};
constexpr riptide::sim::EventMeta kMeta{7, 11, 13};

HadronicStep step(HadronicKind kind, int Z, std::optional<RecoilProton> proton = std::nullopt) {
    HadronicStep s;
    s.kind = kind;
    s.target_Z = Z;
    s.position = {1.0 * u::mm, 2.0 * u::mm, 3.0 * u::mm};
    s.time = 0.5 * u::ns;
    s.proton = proton;
    return s;
}

} // namespace

TEST_CASE("first interaction classification", "[truth]") {
    CHECK(classifyFirstInteraction(HadronicKind::Elastic, kHydrogenZ) ==
          FirstInteraction::NpElastic);
    CHECK(classifyFirstInteraction(HadronicKind::Elastic, kCarbonZ) == FirstInteraction::NCElastic);
    CHECK(classifyFirstInteraction(HadronicKind::Inelastic, kCarbonZ) ==
          FirstInteraction::NCInelastic);
    // Capture on hydrogen must not be mistaken for an n-p scatter.
    CHECK(classifyFirstInteraction(HadronicKind::Capture, kHydrogenZ) == FirstInteraction::Other);
    CHECK(classifyFirstInteraction(HadronicKind::Capture, kCarbonZ) == FirstInteraction::Other);
    CHECK(classifyFirstInteraction(HadronicKind::Inelastic, kHydrogenZ) == FirstInteraction::Other);
    CHECK(classifyFirstInteraction(HadronicKind::Fission, kCarbonZ) == FirstInteraction::Other);
    CHECK(classifyFirstInteraction(HadronicKind::Other, kHydrogenZ) == FirstInteraction::Other);
    CHECK(classifyFirstInteraction(HadronicKind::Elastic, kOxygenZ) == FirstInteraction::Other);
}

TEST_CASE("thetaP on reference vectors", "[truth]") {
    CHECK_THAT(thetaP(kBeam, {0.0, 0.0, 1.0}), WithinAbs(0.0, kAngleTol));
    CHECK_THAT(thetaP(kBeam, {1.0, 0.0, 0.0}), WithinAbs(0.5 * u::pi, kAngleTol));
    CHECK_THAT(thetaP(kBeam, {0.0, 0.0, -1.0}), WithinAbs(u::pi, kAngleTol));
    // Not normalised: only the directions matter.
    CHECK_THAT(thetaP({0.0, 0.0, 3.0}, {2.0, 0.0, 2.0}), WithinRel(0.25 * u::pi, kAngleTol));
}

TEST_CASE("thetaP clamps a cosine rounded above 1", "[truth]") {
    // The normalised dot product of these nearly equal vectors can round to
    // slightly more than 1; the result must be a finite angle near zero.
    const Vec3 a{0.1, 0.2, 0.3};
    const Vec3 b{0.1 * (1.0 + 1e-16), 0.2, 0.3};
    const double theta = thetaP(a, b);
    CHECK(std::isfinite(theta));
    CHECK_THAT(theta, WithinAbs(0.0, 1e-7));
}

TEST_CASE("recorder: new event starts clean", "[truth]") {
    TruthRecorder rec;
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(step(HadronicKind::Elastic, kHydrogenZ, RecoilProton{5.0, {1.0, 0.0, 0.0}}));
    rec.beginEvent({8, 17, 19}, 10.0 * u::MeV, {1.0, 0.0, 0.0});
    CHECK(rec.meta().event_id == 8);
    CHECK(rec.meta().rng_seed_0 == 17);
    CHECK(rec.meta().rng_seed_1 == 19);
    const auto& r = rec.record();
    CHECK_THAT(r.En, WithinAbs(10.0 * u::MeV, 0.0));
    CHECK_THAT(r.omega_n.x, WithinAbs(1.0, 0.0));
    CHECK(r.first_interaction == FirstInteraction::Unknown);
    CHECK(r.n_np_scatters == 0);
    CHECK_FALSE(r.valid_np);
    CHECK(r.n_np_without_proton == 0);
}

TEST_CASE("recorder: n-C elastic then n-p", "[truth]") {
    TruthRecorder rec;
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(step(HadronicKind::Elastic, kCarbonZ));
    rec.onHadronicStep(
        step(HadronicKind::Elastic, kHydrogenZ, RecoilProton{7.0 * u::MeV, {0.6, 0.0, 0.8}}));
    const auto& r = rec.record();
    CHECK(r.first_interaction == FirstInteraction::NCElastic);
    CHECK_THAT(r.r0.z, WithinAbs(3.0 * u::mm, 0.0));
    CHECK_THAT(r.t0, WithinAbs(0.5 * u::ns, 0.0));
    CHECK(r.n_np_scatters == 1);
    REQUIRE(r.valid_np);
    CHECK_THAT(r.Ep, WithinAbs(7.0 * u::MeV, 0.0));
    CHECK_THAT(r.theta_p, WithinRel(std::acos(0.8), kAngleTol));
}

TEST_CASE("recorder: the first n-p defines the proton", "[truth]") {
    TruthRecorder rec;
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(
        step(HadronicKind::Elastic, kHydrogenZ, RecoilProton{7.0 * u::MeV, {0.0, 0.0, 1.0}}));
    rec.onHadronicStep(
        step(HadronicKind::Elastic, kHydrogenZ, RecoilProton{3.0 * u::MeV, {1.0, 0.0, 0.0}}));
    const auto& r = rec.record();
    CHECK(r.first_interaction == FirstInteraction::NpElastic);
    CHECK(r.n_np_scatters == 2);
    CHECK_THAT(r.Ep, WithinAbs(7.0 * u::MeV, 0.0));
    CHECK_THAT(r.theta_p, WithinAbs(0.0, kAngleTol));
}

TEST_CASE("recorder: inelastic gives no n-p", "[truth]") {
    TruthRecorder rec;
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(step(HadronicKind::Inelastic, kCarbonZ));
    const auto& r = rec.record();
    CHECK(r.first_interaction == FirstInteraction::NCInelastic);
    CHECK(r.n_np_scatters == 0);
    CHECK_FALSE(r.valid_np);
}

TEST_CASE("recorder: n-p without a recoil proton", "[truth]") {
    TruthRecorder rec;
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(step(HadronicKind::Elastic, kHydrogenZ));
    // A later n-p with a proton does not replace the first one.
    rec.onHadronicStep(
        step(HadronicKind::Elastic, kHydrogenZ, RecoilProton{3.0 * u::MeV, {1.0, 0.0, 0.0}}));
    const auto& r = rec.record();
    CHECK(r.n_np_scatters == 2);
    CHECK(r.n_np_without_proton == 1);
    CHECK_FALSE(r.valid_np);
}

TEST_CASE("tally counts events per first interaction", "[truth]") {
    TruthRecorder rec;
    TruthTally tally;
    rec.beginEvent(kMeta, kEn, kBeam); // no interaction
    tally.add(rec.record());
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(step(HadronicKind::Elastic, kHydrogenZ)); // n-p, proton missing
    tally.add(rec.record());
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(step(HadronicKind::Inelastic, kCarbonZ));
    tally.add(rec.record());
    CHECK(tally.events() == 3);
    CHECK(tally.count(FirstInteraction::Unknown) == 1);
    CHECK(tally.count(FirstInteraction::NpElastic) == 1);
    CHECK(tally.count(FirstInteraction::NCInelastic) == 1);
    CHECK(tally.count(FirstInteraction::NCElastic) == 0);
    CHECK(tally.count(FirstInteraction::Other) == 0);
    CHECK(tally.npWithoutProton() == 1);
}
