#pragma once

#include <G4UserSteppingAction.hh>

#include "riptide/sim/truth_record.hpp"

class G4Step;

namespace riptide::sim {

// Feeds the truth recorder with the hadronic interactions of the primary
// neutron (track 1) inside the scintillator; every other step is ignored.
class SteppingAction final : public G4UserSteppingAction {
public:
    explicit SteppingAction(TruthRecorder& recorder) : recorder_(&recorder) {}

    void UserSteppingAction(const G4Step* step) override;

private:
    TruthRecorder* recorder_; // owned by the thread's RunAction
};

} // namespace riptide::sim
