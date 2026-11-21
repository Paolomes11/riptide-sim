#pragma once

#include <G4VUserActionInitialization.hh>

#include <cstdint>
#include <filesystem>

#include "riptide/sim/run_config.hpp"

namespace riptide::sim {

class DetectorConstruction;

// Creates the user actions of each worker thread (Build) and of the master
// (BuildForMaster). Each worker gets its own copies: no state is shared.
class ActionInitialization final : public G4VUserActionInitialization {
public:
    // The detector must outlive the run manager's actions; the run manager
    // owns it, so it does.
    ActionInitialization(RunConfig config, std::uint64_t run_seed, std::filesystem::path run_dir,
                         const DetectorConstruction& detector);
    void BuildForMaster() const override;
    void Build() const override;

private:
    RunConfig config_;
    std::uint64_t run_seed_;
    std::filesystem::path run_dir_;
    const DetectorConstruction* detector_;
};

} // namespace riptide::sim
