#include "riptide/sim/action_initialization.hpp"

#include <utility>

#include "riptide/sim/primary_generator_action.hpp"

namespace riptide::sim {

ActionInitialization::ActionInitialization(RunConfig config, std::uint64_t run_seed)
    : config_(std::move(config)), run_seed_(run_seed) {}

void ActionInitialization::BuildForMaster() const {}

void ActionInitialization::Build() const {
    // Ownership passes to the worker run manager.
    SetUserAction(new PrimaryGeneratorAction(config_.source, run_seed_));
}

} // namespace riptide::sim
