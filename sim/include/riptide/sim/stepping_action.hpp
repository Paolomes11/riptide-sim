#pragma once

#include <G4UserSteppingAction.hh>

#include "riptide/sim/deposit_record.hpp"
#include "riptide/sim/truth_record.hpp"

class G4HadronicProcess;
class G4Step;

namespace riptide::sim {

class DetectorConstruction;

// Feeds the recorders with the steps taken inside the scintillator: every
// step goes to the deposit recorder; the hadronic interactions of the primary
// neutron (track 1) and of other neutrons go to the truth recorder. The recoil
// proton of the first n-p scatter is tagged on its track, so its later steps
// are recognised without comparing energies.
class SteppingAction final : public G4UserSteppingAction {
public:
    SteppingAction(TruthRecorder& truth, DepositRecorder& deposits,
                   const DetectorConstruction& detector)
        : truth_(&truth), deposits_(&deposits), detector_(&detector) {}

    void UserSteppingAction(const G4Step* step) override;

private:
    void recordDeposit(const G4Step* step, bool hadronic) const;
    void onPrimaryHadronic(const G4Step* step, const G4HadronicProcess& process) const;

    TruthRecorder* truth_;                 // owned by the thread's RunAction
    DepositRecorder* deposits_;            // owned by the thread's RunAction
    const DetectorConstruction* detector_; // owned by the run manager, shared by all threads
};

} // namespace riptide::sim
