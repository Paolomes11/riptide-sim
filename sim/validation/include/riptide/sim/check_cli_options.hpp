#pragma once

// Command line of riptide_sim_check. Like riptide_sim's, parsing is a pure
// function of the argument list, unit-tested in the CI build; the application
// then resolves the data root and checks the directories.

#include <cstddef>
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

} // namespace riptide::sim::checks
