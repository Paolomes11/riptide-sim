// Entry point of the Geant4 simulation: parses the command line and the run
// configuration, then checks the data root before anything is simulated, so a
// missing SSD or a full disk stops the run at startup rather than mid-run.

#include <G4PhysListFactory.hh>
#include <G4RunManager.hh>
#include <G4RunManagerFactory.hh>
#include <G4VModularPhysicsList.hh>
#include <G4Version.hh>
#include <Randomize.hh>

#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/event_schema.hpp"
#include "riptide/common/version.hpp"
#include "riptide/sim/action_initialization.hpp"
#include "riptide/sim/cli_options.hpp"
#include "riptide/sim/detector_construction.hpp"
#include "riptide/sim/run_config.hpp"

namespace {

// Provisional, deliberately generous estimate of the output size per event
// (uncompressed truth row plus ROOT overhead), used only for the free-space
// check; it is replaced by a measured value in phase 1b.
constexpr std::uint64_t kProvisionalBytesPerEvent = 1024;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;

void printVersion() {
    std::cout << "riptide_sim " << riptide::kProjectVersion << '\n'
              << "Geant4 " << G4VERSION_NUMBER / 100 << '.' << (G4VERSION_NUMBER / 10) % 10 << '.'
              << G4VERSION_NUMBER % 10 << ' ' << G4Version << ' ' << G4Date << '\n'
              << "event schema " << riptide::schema::kSchemaVersion << '\n'
              << "config version " << riptide::sim::kConfigVersion << '\n';
}

// Validates everything that can be validated without Geant4 and returns the run
// directory, which must not exist yet: an existing run is never overwritten.
std::filesystem::path prepareRunDir(const riptide::sim::CliOptions& opts) {
    namespace data = riptide::data;
    const auto root = data::resolveDataRoot(opts.data_dir);
    data::requireWritableRoot(root, data::estimateRunBytes(opts.events, kProvisionalBytesPerEvent),
                              data::kFreeSpaceMargin, opts.allow_system_disk);
    const auto run_dir = data::simRunDir(root, opts.run_id);
    if (std::filesystem::exists(run_dir)) {
        throw data::DataPathError("run directory already exists: " + run_dir.string());
    }
    return run_dir;
}

// Classic MT run manager, requested explicitly: "Default" would let the
// G4RUN_MANAGER_TYPE environment variable switch to another backend.
void runSimulation(const riptide::sim::RunConfig& config, const riptide::sim::CliOptions& opts,
                   const std::filesystem::path& run_dir) {
    G4PhysListFactory factory;
    if (!factory.IsReferencePhysList(config.physics_list)) {
        throw riptide::sim::ConfigError("physics_list: unknown reference physics list '" +
                                        config.physics_list + "'");
    }
    const std::unique_ptr<G4RunManager> run_manager(
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::MTOnly, opts.threads));
    // The per-event seeding relies on MixMax's non-overlapping streams
    // (event_seeds.hpp); workers clone the master's engine type.
    if (dynamic_cast<CLHEP::MixMaxRng*>(G4Random::getTheEngine()) == nullptr) {
        throw std::runtime_error("random engine is not CLHEP::MixMaxRng");
    }
    // The run manager takes ownership of the three initialisation objects.
    run_manager->SetUserInitialization(new riptide::sim::DetectorConstruction(config));
    run_manager->SetUserInitialization(factory.GetReferencePhysList(config.physics_list));
    run_manager->SetUserInitialization(
        new riptide::sim::ActionInitialization(config, opts.seed, run_dir));
    run_manager->Initialize();
    run_manager->BeamOn(static_cast<int>(opts.events));
}

} // namespace

int main(int argc, char** argv) {
    riptide::sim::CliOptions opts;
    try {
        opts = riptide::sim::parseCommandLine(std::vector<std::string>(argv + 1, argv + argc));
    } catch (const riptide::sim::CliError& e) {
        std::cerr << "riptide_sim: " << e.what() << "\n\n" << riptide::sim::helpText();
        return kExitUsage;
    }
    if (opts.show_help) {
        std::cout << riptide::sim::helpText();
        return kExitOk;
    }
    if (opts.show_version) {
        printVersion();
        return kExitOk;
    }

    try {
        const auto config = riptide::sim::loadRunConfig(opts.config);
        const auto run_dir = prepareRunDir(opts);
        std::cout << "riptide_sim: run directory " << run_dir << " (physics list "
                  << config.physics_list << ", " << opts.events << " events, " << opts.threads
                  << " threads, seed " << opts.seed << ")\n";
        // The root was checked by prepareRunDir, so this never writes to the system disk.
        std::filesystem::create_directories(run_dir);
        runSimulation(config, opts, run_dir);
        return kExitOk;
    } catch (const std::exception& e) {
        std::cerr << "riptide_sim: " << e.what() << '\n';
        return kExitRuntime;
    }
}
