#include "riptide/sim/hp_recoil.hpp"

#include <cmath>
#include <string>

#include "riptide/common/np_kinematics.hpp"

namespace riptide::sim::checks {

namespace {

constexpr int kMinRepFlag = 1;
constexpr int kMaxRepFlag = 3;

// sqrt(p^2 + m^2) - m written without the cancellation of the direct form,
// which loses about six digits at a few MeV on a ~1 GeV rest energy.
double kineticFromMomentumSquared(double p_squared, double mass) noexcept {
    return p_squared / (std::sqrt(p_squared + mass * mass) + mass);
}

} // namespace

HpElasticHeader parseHpElasticHeader(std::istream& in) {
    HpElasticHeader h;
    if (!(in >> h.rep_flag >> h.awr >> h.frame_flag)) {
        throw HpHeaderError("HP elastic header: expected 'repFlag targetMass frameFlag'");
    }
    if (h.rep_flag < kMinRepFlag || h.rep_flag > kMaxRepFlag) {
        throw HpHeaderError("HP elastic header: unsupported repFlag " + std::to_string(h.rep_flag));
    }
    if (!std::isfinite(h.awr) || !(h.awr > 0.0)) {
        throw HpHeaderError("HP elastic header: target mass must be positive");
    }
    if (h.frame_flag != kHpFrameCm) {
        throw HpHeaderError("HP elastic header: frameFlag " + std::to_string(h.frame_flag) +
                            ", only the CM frame (2) is modelled");
    }
    return h;
}

double hpProtonEnergy(double En, double theta_p, double target_mass) noexcept {
    const double target_kinetic =
        kinematics::recoilEnergy(En, theta_p, kinematics::kNeutronMass, target_mass);
    const double p_squared = kinematics::momentumSquared(target_kinetic, target_mass);
    return kineticFromMomentumSquared(p_squared, kinematics::kProtonMass);
}

double hpCosCm(double Ep, double En, double target_mass) noexcept {
    const double p_squared = kinematics::momentumSquared(Ep, kinematics::kProtonMass);
    const double target_kinetic = kineticFromMomentumSquared(p_squared, target_mass);
    return 1.0 - 2.0 * target_kinetic /
                     kinematics::maxRecoilEnergy(En, kinematics::kNeutronMass, target_mass);
}

double hpCosCmFromAngle(double En, double theta_p, double target_mass) noexcept {
    const double target_kinetic =
        kinematics::recoilEnergy(En, theta_p, kinematics::kNeutronMass, target_mass);
    return 1.0 - 2.0 * target_kinetic /
                     kinematics::maxRecoilEnergy(En, kinematics::kNeutronMass, target_mass);
}

} // namespace riptide::sim::checks
