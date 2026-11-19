#pragma once

#include <G4UserEventAction.hh>

#include "riptide/sim/deposit_record.hpp"
#include "riptide/sim/ntuple_writer.hpp"
#include "riptide/sim/truth_record.hpp"

class G4Event;

namespace riptide::sim {

// Opens the event's deposit record and closes both records: adds the truth to
// the run tally and writes the row. The truth record is opened by the primary
// generator instead: Geant4 calls GeneratePrimaries before BeginOfEventAction,
// and the generator is the one that knows En and omega_n.
class EventAction final : public G4UserEventAction {
public:
    EventAction(const TruthRecorder& recorder, DepositRecorder& deposits, TruthTally& tally,
                NtupleWriter& writer)
        : recorder_(&recorder), deposits_(&deposits), tally_(&tally), writer_(&writer) {}

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

private:
    const TruthRecorder* recorder_; // owned by the thread's RunAction
    DepositRecorder* deposits_;     // owned by the thread's RunAction
    TruthTally* tally_;             // owned by the thread's RunAction
    NtupleWriter* writer_;          // owned by the thread's RunAction
};

} // namespace riptide::sim
