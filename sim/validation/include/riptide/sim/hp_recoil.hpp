#pragma once

// Header of the G4NDL elastic final-state files, from which the checks read
// the target mass that the Geant4 high-precision model uses (see
// common/hp_kinematics.hpp for why it matters).

#include <istream>
#include <stdexcept>

#include "riptide/common/hp_kinematics.hpp"

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

// The HP kinematics lives in common (hp_kinematics.hpp) so that the analysis
// can remove the convention; the checks keep their unqualified names.
using kinematics::hpCosCm;
using kinematics::hpCosCmFromAngle;
using kinematics::hpProtonEnergy;

} // namespace riptide::sim::checks
