#pragma once

// Relativistic two-body kinematics of elastic scattering on a target at rest,
// specialised to n-p, plus the Legendre polynomials used to compare angular
// distributions with evaluated nuclear data.
//
// All energies are kinetic energies in MeV, angles in rad (riptide::units).
// Formulas are exact relativistic expressions (four-momentum conservation),
// not expansions: at 20 MeV the non-relativistic recoil energy is off by
// ~0.5% at 45 deg, far more than the checks in phase 1b can tolerate.

#include <cmath>

#include "riptide/common/units.hpp"

namespace riptide::kinematics {

// Rest energies m c^2. These are the values of CLHEP (PhysicalConstants.h,
// CODATA 2006) used by Geant4 for G4Proton and G4Neutron, so that recoil
// energies computed here and the ones Geant4 produces share the same
// constants; sim/tests checks the equality. The newer CODATA values differ by
// < 1e-7 relative, below anything the analysis resolves.
inline constexpr double kProtonMass = 938.272013 * units::MeV;
inline constexpr double kNeutronMass = 939.56536 * units::MeV;

// Squared lab momentum (p c)^2 of a projectile of kinetic energy `kinetic` and
// rest energy `mass`.
[[nodiscard]] inline double momentumSquared(double kinetic, double mass) noexcept {
    return kinetic * (kinetic + 2.0 * mass);
}

// The general functions below take (kinetic energy, angle or masses) as plain
// doubles because that is the physics signature; the n-p wrappers further down
// fix the masses and are what the rest of the project calls.
// NOLINTBEGIN(bugprone-easily-swappable-parameters)

// Maximum kinetic energy of the recoiling target (head-on collision, recoil at
// theta = 0): 2 m_t p^2 / s, with s the squared invariant mass. For equal
// masses it reduces to the projectile kinetic energy.
[[nodiscard]] inline double maxRecoilEnergy(double kinetic, double m_projectile,
                                            double m_target) noexcept {
    const double p_squared = momentumSquared(kinetic, m_projectile);
    const double s_invariant = m_projectile * m_projectile + m_target * m_target +
                               2.0 * m_target * (kinetic + m_projectile);
    return 2.0 * m_target * p_squared / s_invariant;
}

// Kinetic energy of the recoiling target at lab angle theta (from the
// projectile direction), theta in [0, pi/2]:
//   T = 2 m_t p^2 cos^2(theta) / ((E + m_t)^2 - p^2 cos^2(theta)),
// with E the projectile total energy. Beyond pi/2 the recoil is kinematically
// forbidden; the formula would still give the mirror value, so callers must
// not rely on it there.
[[nodiscard]] inline double recoilEnergy(double kinetic, double theta, double m_projectile,
                                         double m_target) noexcept {
    const double p_squared = momentumSquared(kinetic, m_projectile);
    const double cos_theta = std::cos(theta);
    const double p2_cos2 = p_squared * cos_theta * cos_theta;
    const double total_plus_target = kinetic + m_projectile + m_target;
    return 2.0 * m_target * p2_cos2 / (total_plus_target * total_plus_target - p2_cos2);
}

// NOLINTEND(bugprone-easily-swappable-parameters)

// Recoil proton energy Ep for a neutron of energy En and proton lab angle theta_p.
[[nodiscard]] inline double protonEnergy(double En, double theta_p) noexcept {
    return recoilEnergy(En, theta_p, kNeutronMass, kProtonMass);
}

// Maximum recoil proton energy for a neutron of energy En.
[[nodiscard]] inline double maxProtonEnergy(double En) noexcept {
    return maxRecoilEnergy(En, kNeutronMass, kProtonMass);
}

// The recoil energy is linear in the CM scattering cosine of the projectile,
// exactly and for any masses: the momentum transfer is t = -2 p*^2 (1 - cos),
// and for a target at rest the recoil energy is -t / (2 m_t). Hence
//   Ep = Ep_max (1 - cos_cm) / 2,
// where cos_cm is the cosine of the neutron's CM angle, the variable mu of
// the evaluated angular distributions (File 4, CM frame).
[[nodiscard]] inline double protonEnergyFromCosCm(double En, double cos_cm) noexcept {
    return 0.5 * maxProtonEnergy(En) * (1.0 - cos_cm);
}

// Inverse of protonEnergyFromCosCm: neutron CM cosine from (Ep, En).
[[nodiscard]] inline double cosCmFromProtonEnergy(double Ep, double En) noexcept {
    return 1.0 - 2.0 * Ep / maxProtonEnergy(En);
}

// Legendre polynomial P_degree(x) by the Bonnet recurrence
//   (k + 1) P_{k+1} = (2k + 1) x P_k - k P_{k-1},
// which is stable for |x| <= 1. A negative degree is not a Legendre degree
// and returns 0.
[[nodiscard]] inline double legendreP(int degree, double x) noexcept {
    if (degree < 0) {
        return 0.0;
    }
    double previous = 1.0; // P_0
    if (degree == 0) {
        return previous;
    }
    double current = x; // P_1
    for (int k = 1; k < degree; ++k) {
        const double next = ((2.0 * k + 1.0) * x * current - k * previous) / (k + 1.0);
        previous = current;
        current = next;
    }
    return current;
}

} // namespace riptide::kinematics
