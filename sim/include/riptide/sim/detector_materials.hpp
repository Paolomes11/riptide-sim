#pragma once

// Names of the Geant4 NIST materials of the detector. The simulation builds
// the scintillator from this name, and riptide_sim_check reads the atom
// densities of the same material, so the expected interaction rates always
// refer to the medium that was actually simulated.

namespace riptide::sim {

// Polyvinyltoluene-based plastic scintillator of the Geant4 NIST database
// (G4NistMaterialBuilder: density 1.032 g/cm3).
inline constexpr const char* kScintillatorMaterial = "G4_PLASTIC_SC_VINYLTOLUENE";

} // namespace riptide::sim
