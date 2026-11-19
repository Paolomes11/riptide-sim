#include "riptide/sim/cli_options.hpp"

#include <lyra/lyra.hpp>

#include <sstream>
#include <string>

#include "riptide/common/data_paths.hpp"
#include "riptide/sim/event_seeds.hpp"

namespace riptide::sim {

namespace {

struct RawOptions {
    bool show_help = false;
    bool show_version = false;
    std::string config;
    std::string data_dir;
    std::string events;
    std::string threads;
    std::string seed;
    std::string run_id;
    bool allow_system_disk = false;
};

lyra::cli makeCli(RawOptions& raw) {
    return lyra::cli() | lyra::help(raw.show_help) |
           lyra::opt(raw.show_version)["--version"]("print the versions and exit") |
           lyra::opt(raw.config, "file")["--config"]("run configuration (JSON), required") |
           lyra::opt(raw.events, "n")["--events"]("number of events, required, > 0") |
           lyra::opt(raw.run_id, "id")["--run-id"]("output directory name under <root>/sim/, "
                                                   "required; must not exist") |
           lyra::opt(raw.data_dir, "dir")["--data-dir"]("data root (default: $RIPTIDE_DATA_DIR, "
                                                        "then /mnt/external_ssd/riptide-sim)") |
           lyra::opt(raw.threads, "n")["--threads"]("worker threads (default 1)") |
           // Built from the constants so the help cannot drift from the check below.
           lyra::opt(raw.seed, "s")["--seed"]("run seed, 0.." + std::to_string(kMaxSeed) +
                                              " (default " + std::to_string(kDefaultSeed) + ")") |
           lyra::opt(raw.allow_system_disk)["--allow-system-disk"](
               "allow writing on the system disk (short tests only)");
}

} // namespace

CliOptions parseCommandLine(const std::vector<std::string>& args) {
    RawOptions raw;
    const auto cli = makeCli(raw);
    std::vector<std::string> argv{"riptide_sim"};
    argv.insert(argv.end(), args.begin(), args.end());
    const auto result = cli.parse(lyra::args(argv.begin(), argv.end()));
    if (!result) {
        throw CliError(result.message());
    }

    CliOptions opts;
    opts.show_help = raw.show_help;
    opts.show_version = raw.show_version;
    if (opts.show_help || opts.show_version) {
        return opts;
    }

    if (raw.config.empty()) {
        throw CliError("--config is required");
    }
    if (raw.events.empty()) {
        throw CliError("--events is required");
    }
    if (raw.run_id.empty()) {
        throw CliError("--run-id is required");
    }
    opts.config = raw.config;
    opts.events = parseCliNumber<std::uint64_t>(raw.events, "--events");
    if (opts.events == 0 || opts.events > kMaxEvents) {
        throw CliError("--events: must be in [1, " + std::to_string(kMaxEvents) + "]");
    }
    if (!data::isValidId(raw.run_id)) {
        throw CliError("--run-id: only [A-Za-z0-9._-], not starting with '.' or '-', at most " +
                       std::to_string(data::kMaxIdLength) + " characters");
    }
    opts.run_id = raw.run_id;
    if (!raw.threads.empty()) {
        opts.threads = parseCliNumber<int>(raw.threads, "--threads");
        if (opts.threads < 1 || opts.threads > kMaxThreads) {
            throw CliError("--threads: must be in [1, " + std::to_string(kMaxThreads) + "]");
        }
    }
    if (!raw.seed.empty()) {
        opts.seed = parseCliNumber<std::uint64_t>(raw.seed, "--seed");
        if (opts.seed > kMaxSeed) {
            throw CliError("--seed: must be in [0, " + std::to_string(kMaxSeed) +
                           "] (seeds are stored as signed 32-bit integers)");
        }
    }
    if (!raw.data_dir.empty()) {
        opts.data_dir = raw.data_dir;
    }
    opts.allow_system_disk = raw.allow_system_disk;
    return opts;
}

std::string helpText() {
    RawOptions raw;
    std::ostringstream out;
    out << makeCli(raw);
    return out.str();
}

} // namespace riptide::sim
