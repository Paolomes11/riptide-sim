#pragma once

#include <G4VUserActionInitialization.hh>

#include "riptide/sim/run_config.hpp"

namespace riptide::sim {

// Creates the user actions of each worker thread (Build) and of the master
// (BuildForMaster). Each worker gets its own copies: no state is shared.
class ActionInitialization final : public G4VUserActionInitialization {
public:
    explicit ActionInitialization(RunConfig config);
    void BuildForMaster() const override;
    void Build() const override;

private:
    RunConfig config_;
};

} // namespace riptide::sim
