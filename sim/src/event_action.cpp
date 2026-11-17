#include "riptide/sim/event_action.hpp"

namespace riptide::sim {

void EventAction::EndOfEventAction(const G4Event* /*event*/) {
    tally_->add(recorder_->record());
}

} // namespace riptide::sim
