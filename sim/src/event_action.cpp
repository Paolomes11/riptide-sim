#include "riptide/sim/event_action.hpp"

namespace riptide::sim {

void EventAction::BeginOfEventAction(const G4Event* /*event*/) {
    // No step of the event has been taken yet: tracking starts after this.
    deposits_->beginEvent();
}

void EventAction::EndOfEventAction(const G4Event* /*event*/) {
    tally_->add(recorder_->record());
    writer_->fill(eventRow(recorder_->meta(), recorder_->record(), deposits_->record()));
}

} // namespace riptide::sim
