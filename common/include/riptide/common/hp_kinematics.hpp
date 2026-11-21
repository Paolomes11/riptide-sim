#pragma once

// Recoil-proton kinematics as the Geant4 high-precision elastic model
// produces it, and the conversion back to the physical n-p kinematics.
//
// G4ParticleHPElasticFS takes the target mass from the final-state file of
// the nuclide: AWR (in neutron masses) on its first line. For H-1 that AWR is
// the mass of the hydrogen atom, electron included, not of the bare proton.
// The scattering is solved with that mass, and the recoil is then created as
// a G4Proton with the resulting momentum, so its kinetic energy is
// sqrt(p^2 + m_p^2) - m_p: higher than the bare-proton kinematics by about
// m_e / m_p in relative terms. Above the thermal (SVT) threshold, a few eV
// at room temperature, the target is at rest, which is what is modelled here.
//
// The simulation checks use these functions to compare Geant4 with what it is
// meant to do; the analysis uses physicalProtonEnergyFromHp to remove the
// convention from the truth, since real protons do not carry it. The target
// mass is always an argument: it is read from the nuclear data the run used,
// never written by hand.

#include <cmath>

#include "riptide/common/np_kinematics.hpp"

namespace riptide::kinematics {

// The functions take (energy, angle or energy, mass) as plain doubles like
// np_kinematics.hpp.
// NOLINTBEGIN(bugprone-easily-swappable-parameters)

// sqrt(p^2 + m^2) - m written without the cancellation of the direct form,
// which loses about six digits at a few MeV on a ~1 GeV rest energy.
[[nodiscard]] inline double kineticFromMomentumSquared(double p_squared, double mass) noexcept {
    return p_squared / (std::sqrt(p_squared + mass * mass) + mass);
}

// Kinetic energy [MeV] of the recoil proton at lab angle theta_p for a
// neutron of energy En on a target of rest energy target_mass at rest.
[[nodiscard]] inline double hpProtonEnergy(double En, double theta_p, double target_mass) noexcept {
    const double target_kinetic = recoilEnergy(En, theta_p, kNeutronMass, target_mass);
    const double p_squared = momentumSquared(target_kinetic, target_mass);
    return kineticFromMomentumSquared(p_squared, kProtonMass);
}

// Neutron CM cosine from the recoil-proton energy: the exact inverse of
// hpProtonEnergy through the target's recoil energy, which is linear in the
// CM cosine for any masses.
[[nodiscard]] inline double hpCosCm(double Ep, double En, double target_mass) noexcept {
    const double p_squared = momentumSquared(Ep, kProtonMass);
    const double target_kinetic = kineticFromMomentumSquared(p_squared, target_mass);
    return 1.0 - 2.0 * target_kinetic / maxRecoilEnergy(En, kNeutronMass, target_mass);
}

// Neutron CM cosine from the proton lab angle with the same target mass.
[[nodiscard]] inline double hpCosCmFromAngle(double En, double theta_p,
                                             double target_mass) noexcept {
    const double target_kinetic = recoilEnergy(En, theta_p, kNeutronMass, target_mass);
    return 1.0 - 2.0 * target_kinetic / maxRecoilEnergy(En, kNeutronMass, target_mass);
}

// Physical recoil-proton energy for the same scattering as a Geant4 HP proton
// of energy Ep_hp: the CM cosine is recovered exactly with the HP target mass
// and mapped back with the bare-proton kinematics. Recomputing from the lab
// angle instead would leave a few 1e-6 of the energy, since the HP model also
// shifts theta_p slightly.
[[nodiscard]] inline double physicalProtonEnergyFromHp(double Ep_hp, double En,
                                                       double target_mass) noexcept {
    return protonEnergyFromCosCm(En, hpCosCm(Ep_hp, En, target_mass));
}

// NOLINTEND(bugprone-easily-swappable-parameters)

} // namespace riptide::kinematics
