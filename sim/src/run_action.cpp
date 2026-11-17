#include "riptide/sim/run_action.hpp"

#include <G4ios.hh>

namespace riptide::sim {

void RunAction::BeginOfRunAction(const G4Run* /*run*/) {
    tally_ = TruthTally{};
}

void RunAction::EndOfRunAction(const G4Run* /*run*/) {
    using schema::FirstInteraction;
    // G4cout on a worker is prefixed with the thread id, so each line is one
    // thread's share of the run.
    G4cout << "truth summary: " << tally_.events() << " events, no interaction "
           << tally_.count(FirstInteraction::Unknown) << ", np_elastic "
           << tally_.count(FirstInteraction::NpElastic) << ", nC_elastic "
           << tally_.count(FirstInteraction::NCElastic) << ", nC_inelastic "
           << tally_.count(FirstInteraction::NCInelastic) << ", other "
           << tally_.count(FirstInteraction::Other) << ", H elastic without recoil proton "
           << tally_.npWithoutProton() << G4endl;
}

} // namespace riptide::sim
