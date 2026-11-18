#pragma once

// Recoil-proton kinematics as the Geant4 high-precision elastic model
// produces it, so that the checks compare the simulation with what it is
// meant to do rather than with a slightly different convention.
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
// common/np_kinematics.hpp stays the physical n-p kinematics used by the
// analysis; these functions describe the simulation only.

#include <istream>
#include <stdexcept>

namespace riptide::sim::checks {

class HpHeaderError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// First line of a G4NDL elastic final-state file: representation flag,
// target mass in neutron masses and frame of the angular data.
struct HpElasticHeader {
    int rep_flag = 0;
    double awr = 0.0;
    int frame_flag = 0;
};

// The angular data are tabulated in the CM frame; the kinematics below
// assume it, and a file in the lab frame (1) is refused.
inline constexpr int kHpFrameCm = 2;

// Reads the header. Throws HpHeaderError on malformed text, on a
// representation flag outside 1-3 (0, isotropic, stores the frame elsewhere),
// on AWR <= 0 or on a frame other than kHpFrameCm.
[[nodiscard]] HpElasticHeader parseHpElasticHeader(std::istream& in);

// The functions take (energy, angle or energy, mass) as plain doubles like
// common/np_kinematics.hpp.
// NOLINTBEGIN(bugprone-easily-swappable-parameters)

// Kinetic energy [MeV] of the recoil proton at lab angle theta_p for a
// neutron of energy En on a target of rest energy target_mass at rest.
[[nodiscard]] double hpProtonEnergy(double En, double theta_p, double target_mass) noexcept;

// Neutron CM cosine from the recoil-proton energy: the exact inverse of
// hpProtonEnergy through the target's recoil energy, which is linear in the
// CM cosine for any masses.
[[nodiscard]] double hpCosCm(double Ep, double En, double target_mass) noexcept;

// Neutron CM cosine from the proton lab angle with the same target mass.
[[nodiscard]] double hpCosCmFromAngle(double En, double theta_p, double target_mass) noexcept;

// NOLINTEND(bugprone-easily-swappable-parameters)

} // namespace riptide::sim::checks
