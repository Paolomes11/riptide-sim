#include "riptide/sim/action_initialization.hpp"

#include <utility>

#include "riptide/sim/primary_generator_action.hpp"

namespace riptide::sim {

ActionInitialization::ActionInitialization(RunConfig config) : config_(std::move(config)) {}

void ActionInitialization::BuildForMaster() const {}

void ActionInitialization::Build() const {
    // Ownership passes to the worker run manager.
    SetUserAction(new PrimaryGeneratorAction(config_.source));
}

} // namespace riptide::sim
