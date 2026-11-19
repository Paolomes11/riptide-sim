#pragma once

// Command lines of riptide_sim_check and riptide_sim_predict. Like
// riptide_sim's, parsing is a pure function of the argument list, unit-tested
// in the CI build; the applications then resolve the data root and check the
// directories.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/cli_parse.hpp"

namespace riptide::sim::checks {

// Upper bound on --bins-u and --bins-z: a guard against typos, not a physics
// limit. Bins far beyond the event count only end up merged again.
inline constexpr std::size_t kMaxBins = 1000;

struct CheckCliOptions {
    bool show_help = false;
    bool show_version = false;
    std::optional<std::string> data_dir; // unset: RIPTIDE_DATA_DIR, then the default root
    std::string run_id;
    std::string check_id;
    std::size_t bins_u = kDefaultBins;
    std::size_t bins_z = kDefaultBins;
};

// `args` excludes the program name. Throws CliError on unknown options,
// malformed numbers, invalid ids and missing required options. With --help or
// --version the other options are not required.
[[nodiscard]] CheckCliOptions parseCheckCommandLine(const std::vector<std::string>& args);

// Usage text printed by --help.
[[nodiscard]] std::string checkHelpText();

// Default number of midpoint nodes for the spectrum average of a flat
// source. The evaluated cross sections of carbon are tabulated with fine
// structure, so a coarse rule misses it: over 10-20 MeV the first-interaction
// fractions move by about 0.9 statistical sigma of a 10^7-event run between 64
// and 128 nodes, and by a few 1e-3 sigma between 4096 and 8192. The prediction
// reports that shift at twice the nodes, so the default must keep it far below
// the errors; the cost is well under a minute.
inline constexpr std::size_t kDefaultEnergyNodes = 4096;
// Guard against typos, like kMaxBins.
inline constexpr std::size_t kMaxEnergyNodes = 100000;

struct PredictCliOptions {
    bool show_help = false;
    bool show_version = false;
    std::optional<std::string> data_dir; // unset: RIPTIDE_DATA_DIR, then the default root
    std::string config_path;
    std::uint64_t events = 0;
    std::size_t bins_u = kDefaultBins;
    std::size_t bins_z = kDefaultBins;
    std::size_t energy_nodes = kDefaultEnergyNodes;
};

// `args` excludes the program name. Throws CliError on unknown options,
// malformed numbers, zero events and missing required options. With --help
// or --version the other options are not required.
[[nodiscard]] PredictCliOptions parsePredictCommandLine(const std::vector<std::string>& args);

[[nodiscard]] std::string predictHelpText();

} // namespace riptide::sim::checks
