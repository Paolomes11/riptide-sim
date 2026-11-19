#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <optional>
#include <stdexcept>

#include "riptide/common/units.hpp"
#include "riptide/sim/truth_record.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::schema::EventClass;
using riptide::schema::FirstInteraction;
using riptide::sim::classifyEvent;
using riptide::sim::classifyFirstInteraction;
using riptide::sim::HadronicKind;
using riptide::sim::HadronicStep;
using riptide::sim::InteractionCounts;
using riptide::sim::interactionCounts;
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

TEST_CASE("event classification", "[truth]") {
    CHECK(classifyEvent({0, 0, 0, 0}) == EventClass::NoInteraction);
    CHECK(classifyEvent({1, 0, 0, 0}) == EventClass::SingleNp);
    CHECK(classifyEvent({2, 0, 0, 0}) == EventClass::MultipleNp);
    CHECK(classifyEvent({5, 0, 0, 0}) == EventClass::MultipleNp);
    CHECK(classifyEvent({0, 1, 0, 0}) == EventClass::NCElasticOnly);
    CHECK(classifyEvent({0, 3, 0, 0}) == EventClass::NCElasticOnly);
    CHECK(classifyEvent({0, 0, 1, 0}) == EventClass::NCInelasticNoNp);
    CHECK(classifyEvent({0, 2, 1, 0}) == EventClass::NCInelasticNoNp);
    CHECK(classifyEvent({1, 1, 0, 0}) == EventClass::Mixed);
    CHECK(classifyEvent({1, 0, 1, 0}) == EventClass::Mixed);
    CHECK(classifyEvent({3, 2, 1, 0}) == EventClass::Mixed);
    // Any other interaction wins over everything else.
    CHECK(classifyEvent({0, 0, 0, 1}) == EventClass::Other);
    CHECK(classifyEvent({1, 0, 0, 1}) == EventClass::Other);
    CHECK(classifyEvent({1, 1, 1, 1}) == EventClass::Other);
}

TEST_CASE("event classification refuses negative counts", "[truth]") {
    CHECK_THROWS_AS(classifyEvent({-1, 0, 0, 0}), std::invalid_argument);
    CHECK_THROWS_AS(classifyEvent({0, -1, 0, 0}), std::invalid_argument);
    CHECK_THROWS_AS(classifyEvent({0, 0, -1, 0}), std::invalid_argument);
    CHECK_THROWS_AS(classifyEvent({0, 0, 0, -1}), std::invalid_argument);
}

TEST_CASE("recorder counts every kind of primary interaction", "[truth]") {
    TruthRecorder rec;
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(step(HadronicKind::Elastic, kCarbonZ));
    rec.onHadronicStep(step(HadronicKind::Elastic, kCarbonZ));
    rec.onHadronicStep(
        step(HadronicKind::Elastic, kHydrogenZ, RecoilProton{7.0 * u::MeV, {0.6, 0.0, 0.8}}));
    rec.onHadronicStep(step(HadronicKind::Capture, kHydrogenZ));
    rec.onHadronicStep(step(HadronicKind::Inelastic, kCarbonZ));
    const auto& r = rec.record();
    CHECK(r.n_np_scatters == 1);
    CHECK(r.n_nC_elastic == 2);
    CHECK(r.n_nC_inelastic == 1);
    CHECK(r.n_other == 1);
    const InteractionCounts c = interactionCounts(r);
    CHECK(c.np == 1);
    CHECK(c.nC_elastic == 2);
    CHECK(c.nC_inelastic == 1);
    CHECK(c.other == 1);
    CHECK(classifyEvent(c) == EventClass::Other);
}

TEST_CASE("recorder: n-p vertex is the first n-p, not r0", "[truth]") {
    TruthRecorder rec;
    rec.beginEvent(kMeta, kEn, kBeam);
    auto carbon = step(HadronicKind::Elastic, kCarbonZ);
    carbon.position = {-4.0 * u::mm, 0.0, 0.0};
    rec.onHadronicStep(carbon);
    auto first_np =
        step(HadronicKind::Elastic, kHydrogenZ, RecoilProton{7.0 * u::MeV, {0.0, 0.0, 1.0}});
    first_np.position = {9.0 * u::mm, 8.0 * u::mm, 7.0 * u::mm};
    rec.onHadronicStep(first_np);
    rec.onHadronicStep(
        step(HadronicKind::Elastic, kHydrogenZ, RecoilProton{2.0 * u::MeV, {0.0, 1.0, 0.0}}));
    const auto& r = rec.record();
    CHECK_THAT(r.r0.x, WithinAbs(-4.0 * u::mm, 0.0));
    CHECK_THAT(r.np_vertex.x, WithinAbs(9.0 * u::mm, 0.0));
    CHECK_THAT(r.np_vertex.y, WithinAbs(8.0 * u::mm, 0.0));
    CHECK_THAT(r.np_vertex.z, WithinAbs(7.0 * u::mm, 0.0));
    CHECK(classifyEvent(interactionCounts(r)) == EventClass::Mixed);
}

TEST_CASE("recorder: secondary-neutron interactions are counted apart", "[truth]") {
    TruthRecorder rec;
    rec.beginEvent(kMeta, kEn, kBeam);
    rec.onHadronicStep(step(HadronicKind::Inelastic, kCarbonZ));
    rec.onSecondaryNeutronHadronic();
    rec.onSecondaryNeutronHadronic();
    const auto& r = rec.record();
    CHECK(r.n_secondary_n_hadronic == 2);
    // They do not change the class, which is about the primary only.
    CHECK(classifyEvent(interactionCounts(r)) == EventClass::NCInelasticNoNp);
    rec.beginEvent(kMeta, kEn, kBeam);
    CHECK(rec.record().n_secondary_n_hadronic == 0);
    CHECK(rec.record().n_nC_inelastic == 0);
}
