#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sstream>
#include <string>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/hp_recoil.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace riptide::sim::checks;
namespace kin = riptide::kinematics;
namespace units = riptide::units;

// Target mass of H-1 in the G4NDL 4.7.1 elastic file, in neutron masses.
constexpr double kAwrH1 = 0.999167;
constexpr double kTargetMass = kAwrH1 * kin::kNeutronMass;

constexpr double kEnMin = 10.0 * units::MeV;
constexpr double kEnMax = 20.0 * units::MeV;
constexpr double kEn14 = 14.0 * units::MeV;
constexpr int kGridSteps = 40;

// Equal masses agree to the last few ulps; the round trip loses a few more
// through the square roots; head-on is exact up to rounding.
constexpr double kSameMassRel = 1.0e-14;
constexpr double kRoundTripAbs = 1.0e-12;
constexpr double kHeadOnAbs = 1.0e-12;

// Relative excess of the HP proton energy over the bare-proton kinematics at
// 14 MeV over the open angular range, computed independently from the same
// masses.
constexpr double kExcessMin14 = 5.36e-4;
constexpr double kExcessMax14 = 5.49e-4;

double gridAngle(int i) {
    return 0.5 * units::pi * static_cast<double>(i) / static_cast<double>(kGridSteps + 1);
}

} // namespace

TEST_CASE("with the bare proton mass the HP kinematics is the physical one", "[hp_recoil]") {
    for (int i = 1; i <= kGridSteps; ++i) {
        const double theta = gridAngle(i);
        CHECK_THAT(hpProtonEnergy(kEn14, theta, kin::kProtonMass),
                   WithinRel(kin::protonEnergy(kEn14, theta), kSameMassRel));
    }
}

TEST_CASE("hpCosCm inverts hpProtonEnergy and matches the angle map", "[hp_recoil]") {
    for (int j = 0; j <= kGridSteps; ++j) {
        const double En =
            kEnMin + (kEnMax - kEnMin) * static_cast<double>(j) / static_cast<double>(kGridSteps);
        for (int i = 1; i <= kGridSteps; ++i) {
            const double theta = gridAngle(i);
            const double Ep = hpProtonEnergy(En, theta, kTargetMass);
            CHECK_THAT(hpCosCm(Ep, En, kTargetMass),
                       WithinAbs(hpCosCmFromAngle(En, theta, kTargetMass), kRoundTripAbs));
        }
    }
}

TEST_CASE("a head-on recoil maps to cos_cm = -1", "[hp_recoil]") {
    const double Ep = hpProtonEnergy(kEn14, 0.0, kTargetMass);
    CHECK_THAT(hpCosCm(Ep, kEn14, kTargetMass), WithinAbs(-1.0, kHeadOnAbs));
    CHECK_THAT(hpCosCmFromAngle(kEn14, 0.0, kTargetMass), WithinAbs(-1.0, kHeadOnAbs));
}

TEST_CASE("the atomic target mass raises Ep by about m_e / m_p", "[hp_recoil]") {
    for (int i = 1; i <= kGridSteps; ++i) {
        const double theta = gridAngle(i);
        const double excess =
            hpProtonEnergy(kEn14, theta, kTargetMass) / kin::protonEnergy(kEn14, theta) - 1.0;
        CHECK(excess > kExcessMin14);
        CHECK(excess < kExcessMax14);
    }
}

TEST_CASE("parseHpElasticHeader reads the G4NDL first line", "[hp_recoil]") {
    std::istringstream in(" \t1 \t0.999167 \t2\n \t96\n");
    const HpElasticHeader h = parseHpElasticHeader(in);
    CHECK(h.rep_flag == 1);
    CHECK_THAT(h.awr, WithinRel(kAwrH1, kSameMassRel));
    CHECK(h.frame_flag == kHpFrameCm);
}

TEST_CASE("parseHpElasticHeader refuses what the model does not cover", "[hp_recoil]") {
    for (const std::string text : {"1 0.999167 1", "1 -0.5 2", "0 0.999167 2", "1 0.999167", ""}) {
        std::istringstream in(text);
        CHECK_THROWS_AS(parseHpElasticHeader(in), HpHeaderError);
    }
}
