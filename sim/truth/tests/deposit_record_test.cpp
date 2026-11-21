#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "riptide/common/units.hpp"
#include "riptide/sim/deposit_record.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::sim::DepositRecorder;
using riptide::sim::DepositStep;
using riptide::sim::TrackRole;
using riptide::sim::Vec3;
namespace u = riptide::units;

constexpr double kRelTol = 1e-12;

DepositStep step(TrackRole role, double edep, double evis, const Vec3& pre, const Vec3& post,
                 double length) {
    DepositStep s;
    s.role = role;
    s.edep = edep;
    s.evis = evis;
    s.pre = pre;
    s.post = post;
    s.length = length;
    return s;
}

} // namespace

TEST_CASE("deposit: empty event", "[deposit]") {
    DepositRecorder rec;
    rec.beginEvent();
    const auto& r = rec.record();
    CHECK_THAT(r.Edep, WithinAbs(0.0, 0.0));
    CHECK_THAT(r.Evis, WithinAbs(0.0, 0.0));
    CHECK_FALSE(r.proton.seen);
    CHECK_FALSE(r.proton.exits);
    CHECK(r.proton.n_hadronic == 0);
    CHECK(r.proton.n_nuclear_recoil == 0);
}

TEST_CASE("deposit: totals sum every track, the proton only its own", "[deposit]") {
    DepositRecorder rec;
    rec.beginEvent();
    rec.onStep(step(TrackRole::RecoilProton, 2.0 * u::MeV, 1.5 * u::MeV, {0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.3 * u::mm}, 0.3 * u::mm));
    rec.onStep(step(TrackRole::Other, 0.25 * u::MeV, 0.2 * u::MeV, {5.0, 5.0, 5.0}, {5.0, 5.0, 6.0},
                    1.0 * u::mm));
    rec.onStep(step(TrackRole::RecoilProton, 3.0 * u::MeV, 1.0 * u::MeV, {0.0, 0.0, 0.3 * u::mm},
                    {0.0, 0.1 * u::mm, 0.5 * u::mm}, 0.2 * u::mm));
    const auto& r = rec.record();
    CHECK_THAT(r.Edep, WithinRel(5.25 * u::MeV, kRelTol));
    CHECK_THAT(r.Evis, WithinRel(2.7 * u::MeV, kRelTol));
    REQUIRE(r.proton.seen);
    CHECK_THAT(r.proton.Edep, WithinRel(5.0 * u::MeV, kRelTol));
    CHECK_THAT(r.proton.Evis, WithinRel(2.5 * u::MeV, kRelTol));
    CHECK_THAT(r.proton.length, WithinRel(0.5 * u::mm, kRelTol));
    // Start is the first pre-step point, end the last post-step point.
    CHECK_THAT(r.proton.start.z, WithinAbs(0.0, 0.0));
    CHECK_THAT(r.proton.end.y, WithinAbs(0.1 * u::mm, 0.0));
    CHECK_THAT(r.proton.end.z, WithinAbs(0.5 * u::mm, 0.0));
    CHECK_FALSE(r.proton.exits);
}

TEST_CASE("deposit: proton leaving the volume", "[deposit]") {
    DepositRecorder rec;
    rec.beginEvent();
    auto s = step(TrackRole::RecoilProton, 1.0 * u::MeV, 0.5 * u::MeV, {0.0, 0.0, 49.0 * u::mm},
                  {0.0, 0.0, 50.0 * u::mm}, 1.0 * u::mm);
    s.leaves_volume = true;
    rec.onStep(s);
    CHECK(rec.record().proton.exits);
    CHECK_THAT(rec.record().proton.end.z, WithinAbs(50.0 * u::mm, 0.0));
}

TEST_CASE("deposit: hadronic interactions of the proton are counted", "[deposit]") {
    DepositRecorder rec;
    rec.beginEvent();
    auto s = step(TrackRole::RecoilProton, 1.0 * u::MeV, 0.5 * u::MeV, {}, {}, 1.0 * u::mm);
    s.hadronic = true;
    rec.onStep(s);
    // A hadronic step of another track is not the proton's.
    s.role = TrackRole::Other;
    rec.onStep(s);
    CHECK(rec.record().proton.n_hadronic == 1);
}

TEST_CASE("deposit: nuclear recoils of the proton are counted", "[deposit]") {
    DepositRecorder rec;
    rec.beginEvent();
    auto s = step(TrackRole::RecoilProton, 1.0 * u::MeV, 0.5 * u::MeV, {}, {}, 1.0 * u::mm);
    s.nuclear_recoil = true;
    rec.onStep(s);
    rec.onStep(s);
    s.role = TrackRole::Other;
    rec.onStep(s);
    CHECK(rec.record().proton.n_nuclear_recoil == 2);
    CHECK(rec.record().proton.n_hadronic == 0);
    rec.beginEvent();
    CHECK(rec.record().proton.n_nuclear_recoil == 0);
}

TEST_CASE("deposit: new event starts clean", "[deposit]") {
    DepositRecorder rec;
    rec.beginEvent();
    auto s = step(TrackRole::RecoilProton, 1.0 * u::MeV, 0.5 * u::MeV, {}, {}, 1.0 * u::mm);
    s.leaves_volume = true;
    s.hadronic = true;
    rec.onStep(s);
    rec.beginEvent();
    const auto& r = rec.record();
    CHECK_THAT(r.Edep, WithinAbs(0.0, 0.0));
    CHECK_FALSE(r.proton.seen);
    CHECK_FALSE(r.proton.exits);
    CHECK(r.proton.n_hadronic == 0);
    CHECK_THAT(r.proton.length, WithinAbs(0.0, 0.0));
}
