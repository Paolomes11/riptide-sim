#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "riptide/common/hp_kinematics.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
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

double gridEnergy(int j) {
    return kEnMin + (kEnMax - kEnMin) * static_cast<double>(j) / static_cast<double>(kGridSteps);
}

} // namespace

TEST_CASE("with the bare proton mass the HP kinematics is the physical one", "[hp_kinematics]") {
    for (int i = 1; i <= kGridSteps; ++i) {
        const double theta = gridAngle(i);
        CHECK_THAT(kin::hpProtonEnergy(kEn14, theta, kin::kProtonMass),
                   WithinRel(kin::protonEnergy(kEn14, theta), kSameMassRel));
    }
}

TEST_CASE("hpCosCm inverts hpProtonEnergy and matches the angle map", "[hp_kinematics]") {
    for (int j = 0; j <= kGridSteps; ++j) {
        const double En = gridEnergy(j);
        for (int i = 1; i <= kGridSteps; ++i) {
            const double theta = gridAngle(i);
            const double Ep = kin::hpProtonEnergy(En, theta, kTargetMass);
            CHECK_THAT(kin::hpCosCm(Ep, En, kTargetMass),
                       WithinAbs(kin::hpCosCmFromAngle(En, theta, kTargetMass), kRoundTripAbs));
        }
    }
}

TEST_CASE("a head-on recoil maps to cos_cm = -1", "[hp_kinematics]") {
    const double Ep = kin::hpProtonEnergy(kEn14, 0.0, kTargetMass);
    CHECK_THAT(kin::hpCosCm(Ep, kEn14, kTargetMass), WithinAbs(-1.0, kHeadOnAbs));
    CHECK_THAT(kin::hpCosCmFromAngle(kEn14, 0.0, kTargetMass), WithinAbs(-1.0, kHeadOnAbs));
}

TEST_CASE("the atomic target mass raises Ep by about m_e / m_p", "[hp_kinematics]") {
    for (int i = 1; i <= kGridSteps; ++i) {
        const double theta = gridAngle(i);
        const double excess =
            kin::hpProtonEnergy(kEn14, theta, kTargetMass) / kin::protonEnergy(kEn14, theta) - 1.0;
        CHECK(excess > kExcessMin14);
        CHECK(excess < kExcessMax14);
    }
}
