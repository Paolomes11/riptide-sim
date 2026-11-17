#pragma once

// Command line of riptide_sim. Parsing is a pure function of the argument list
// (no Geant4, no filesystem access), so it is unit-tested in the CI build; the
// application then validates paths and free space with riptide::data.

#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace riptide::sim {

// Common seed of the project, used when --seed is not given.
inline constexpr std::uint64_t kDefaultSeed = 20260907;
// Upper bound on --threads: a guard against typos (e.g. 160 instead of 16), not
// a statement about the machine; the application also checks the core count.
inline constexpr int kMaxThreads = 256;

class CliError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct CliOptions {
    bool show_help = false;
    bool show_version = false;
    std::filesystem::path config;
    std::optional<std::string> data_dir; // unset: RIPTIDE_DATA_DIR, then the default root
    std::uint64_t events = 0;
    int threads = 1;
    std::uint64_t seed = kDefaultSeed;
    std::string run_id;
    // Only for short tests (decision D5): lets the run write on the system disk.
    bool allow_system_disk = false;
};

// `args` excludes the program name. Throws CliError on unknown options,
// malformed numbers and missing required options. With --help or --version
// the other options are not required.
[[nodiscard]] CliOptions parseCommandLine(const std::vector<std::string>& args);

// Usage text printed by --help.
[[nodiscard]] std::string helpText();

} // namespace riptide::sim
