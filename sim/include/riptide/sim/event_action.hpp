#pragma once

#include <G4UserEventAction.hh>

#include "riptide/sim/ntuple_writer.hpp"
#include "riptide/sim/truth_record.hpp"

class G4Event;

namespace riptide::sim {

// Closes the event's truth record: adds it to the run tally and writes its row. The record is
// opened by the primary generator, not here: Geant4 calls GeneratePrimaries before
// BeginOfEventAction.
class EventAction final : public G4UserEventAction {
public:
    EventAction(const TruthRecorder& recorder, TruthTally& tally, NtupleWriter& writer)
        : recorder_(&recorder), tally_(&tally), writer_(&writer) {}

    void EndOfEventAction(const G4Event* event) override;

private:
    const TruthRecorder* recorder_; // owned by the thread's RunAction
    TruthTally* tally_;             // owned by the thread's RunAction
    NtupleWriter* writer_;          // owned by the thread's RunAction
};

} // namespace riptide::sim
