#pragma once

#include <G4VUserActionInitialization.hh>

#include <cstdint>

#include "riptide/sim/run_config.hpp"

namespace riptide::sim {

// Creates the user actions of each worker thread (Build) and of the master
// (BuildForMaster). Each worker gets its own copies: no state is shared.
class ActionInitialization final : public G4VUserActionInitialization {
public:
    ActionInitialization(RunConfig config, std::uint64_t run_seed);
    void BuildForMaster() const override;
    void Build() const override;

private:
    RunConfig config_;
    std::uint64_t run_seed_;
};

} // namespace riptide::sim
