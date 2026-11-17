// Entry point of the Geant4 simulation: parses the command line and the run
// configuration, then checks the data root before anything is simulated, so a
// missing SSD or a full disk stops the run at startup rather than mid-run.

#include <G4Version.hh>

#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/event_schema.hpp"
#include "riptide/common/version.hpp"
#include "riptide/sim/cli_options.hpp"
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
        // The Geant4 run is added in the next step; until then the application
        // stops after the checks without creating anything.
        std::cout << "riptide_sim: configuration and data root OK, run directory " << run_dir
                  << " (physics list " << config.physics_list << ", " << opts.events << " events, "
                  << opts.threads << " threads, seed " << opts.seed << ")\n"
                  << "riptide_sim: the Geant4 run is not implemented yet\n";
        return kExitRuntime;
    } catch (const std::exception& e) {
        std::cerr << "riptide_sim: " << e.what() << '\n';
        return kExitRuntime;
    }
}
