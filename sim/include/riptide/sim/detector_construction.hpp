#pragma once

#include <G4VUserDetectorConstruction.hh>

#include "riptide/sim/run_config.hpp"

class G4VPhysicalVolume;

namespace riptide::sim {

// Volume names, also used by the truth recording to recognise the scintillator.
inline constexpr const char* kWorldName = "World";
inline constexpr const char* kScintillatorName = "Scintillator";

// A single plastic scintillator box centred at the origin, in a vacuum world:
// with no material around it the neutron cannot interact before reaching the
// box, so the first interaction is always inside the detector or nowhere.
class DetectorConstruction final : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(RunConfig config);
    G4VPhysicalVolume* Construct() override;

private:
    RunConfig config_;
};

} // namespace riptide::sim
