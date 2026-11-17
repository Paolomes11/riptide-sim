#pragma once

#include <G4UserRunAction.hh>

#include "riptide/sim/truth_record.hpp"

class G4Run;

namespace riptide::sim {

// Worker-thread run action. It owns the thread's truth recorder (shared by
// reference with the generator, stepping and event actions) and the run
// summary printed at the end of the run.
class RunAction final : public G4UserRunAction {
public:
    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

    [[nodiscard]] TruthRecorder& recorder() {
        return recorder_;
    }
    [[nodiscard]] TruthTally& tally() {
        return tally_;
    }

private:
    TruthRecorder recorder_;
    TruthTally tally_;
};

} // namespace riptide::sim
