#include "riptide/sim/action_initialization.hpp"

#include <utility>

#include "riptide/sim/detector_construction.hpp"
#include "riptide/sim/event_action.hpp"
#include "riptide/sim/primary_generator_action.hpp"
#include "riptide/sim/run_action.hpp"
#include "riptide/sim/stepping_action.hpp"

namespace riptide::sim {

ActionInitialization::ActionInitialization(RunConfig config, std::uint64_t run_seed,
                                           std::filesystem::path run_dir,
                                           const DetectorConstruction& detector)
    : config_(std::move(config)), run_seed_(run_seed), run_dir_(std::move(run_dir)),
      detector_(&detector) {}

// The master gets no actions. In particular it never creates an analysis
// manager, so it writes no file: only the workers write their ntuples.
void ActionInitialization::BuildForMaster() const {}

void ActionInitialization::Build() const {
    // Ownership passes to the worker run manager, which outlives every event, so
    // the other actions may keep references into the RunAction's truth state.
    auto* run_action = new RunAction(run_dir_);
    SetUserAction(run_action);
    SetUserAction(new EventAction(run_action->recorder(), run_action->deposits(),
                                  run_action->tally(), run_action->writer()));
    SetUserAction(new SteppingAction(run_action->recorder(), run_action->deposits(), *detector_));
    SetUserAction(new PrimaryGeneratorAction(config_.source, run_seed_, run_action->recorder()));
}

} // namespace riptide::sim
