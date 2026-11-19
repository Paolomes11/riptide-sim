#include "riptide/sim/check_cli_options.hpp"

#include <lyra/lyra.hpp>

#include <sstream>
#include <string>
#include <string_view>

#include "riptide/common/data_paths.hpp"

namespace riptide::sim::checks {

namespace {

struct RawOptions {
    bool show_help = false;
    bool show_version = false;
    std::string data_dir;
    std::string run_id;
    std::string check_id;
    std::string bins_u;
    std::string bins_z;
};

lyra::cli makeCli(RawOptions& raw) {
    // Built from the constants so the help cannot drift from the checks below.
    const std::string bins_help = "number of bins, 1.." + std::to_string(kMaxBins) + " (default " +
                                  std::to_string(kDefaultBins) + ")";
    return lyra::cli() | lyra::help(raw.show_help) |
           lyra::opt(raw.show_version)["--version"]("print the versions and exit") |
           lyra::opt(raw.run_id, "id")["--run-id"]("run to check, under <root>/sim/, required") |
           lyra::opt(raw.check_id, "id")["--check-id"](
               "output directory name under <root>/validation/<run-id>/, required; must not "
               "exist") |
           lyra::opt(raw.data_dir, "dir")["--data-dir"]("data root (default: $RIPTIDE_DATA_DIR, "
                                                        "then /mnt/external_ssd/riptide-sim)") |
           lyra::opt(raw.bins_u, "n")["--bins-u"]("u = Ep / Ep_max histogram: " + bins_help) |
           lyra::opt(raw.bins_z, "n")["--bins-z"]("depth histogram: " + bins_help);
}

std::string requireId(const std::string& value, std::string_view option) {
    if (value.empty()) {
        throw CliError(std::string{option} + " is required");
    }
    if (!data::isValidId(value)) {
        throw CliError(std::string{option} +
                       ": only [A-Za-z0-9._-], not starting with '.' or '-', at most " +
                       std::to_string(data::kMaxIdLength) + " characters");
    }
    return value;
}

struct RawPredictOptions {
    bool show_help = false;
    bool show_version = false;
    std::string data_dir;
    std::string config_path;
    std::string events;
    std::string bins_u;
    std::string bins_z;
    std::string energy_nodes;
};

lyra::cli makePredictCli(RawPredictOptions& raw) {
    const std::string bins_help = "number of bins, 1.." + std::to_string(kMaxBins) + " (default " +
                                  std::to_string(kDefaultBins) + ")";
    return lyra::cli() | lyra::help(raw.show_help) |
           lyra::opt(raw.show_version)["--version"]("print the versions and exit") |
           lyra::opt(raw.config_path, "path")["--config"]("run configuration, required") |
           lyra::opt(raw.events, "n")["--events"]("events of the planned run, required") |
           lyra::opt(raw.data_dir, "dir")["--data-dir"]("data root (default: $RIPTIDE_DATA_DIR, "
                                                        "then /mnt/external_ssd/riptide-sim)") |
           lyra::opt(raw.bins_u, "n")["--bins-u"]("u = Ep / Ep_max histogram: " + bins_help) |
           lyra::opt(raw.bins_z, "n")["--bins-z"]("depth histogram: " + bins_help) |
           lyra::opt(raw.energy_nodes, "n")["--nodes"](
               "midpoint nodes of the spectrum average, 1.." + std::to_string(kMaxEnergyNodes) +
               " (default " + std::to_string(kDefaultEnergyNodes) + ")");
}

std::size_t parseBins(const std::string& value, std::string_view option) {
    if (value.empty()) {
        return kDefaultBins;
    }
    const auto bins = parseCliNumber<std::size_t>(value, option);
    if (bins == 0 || bins > kMaxBins) {
        throw CliError(std::string{option} + ": must be in [1, " + std::to_string(kMaxBins) + "]");
    }
    return bins;
}

} // namespace

CheckCliOptions parseCheckCommandLine(const std::vector<std::string>& args) {
    RawOptions raw;
    const auto cli = makeCli(raw);
    std::vector<std::string> argv{"riptide_sim_check"};
    argv.insert(argv.end(), args.begin(), args.end());
    const auto result = cli.parse(lyra::args(argv.begin(), argv.end()));
    if (!result) {
        throw CliError(result.message());
    }

    CheckCliOptions opts;
    opts.show_help = raw.show_help;
    opts.show_version = raw.show_version;
    if (opts.show_help || opts.show_version) {
        return opts;
    }
    opts.run_id = requireId(raw.run_id, "--run-id");
    opts.check_id = requireId(raw.check_id, "--check-id");
    opts.bins_u = parseBins(raw.bins_u, "--bins-u");
    opts.bins_z = parseBins(raw.bins_z, "--bins-z");
    if (!raw.data_dir.empty()) {
        opts.data_dir = raw.data_dir;
    }
    return opts;
}

PredictCliOptions parsePredictCommandLine(const std::vector<std::string>& args) {
    RawPredictOptions raw;
    const auto cli = makePredictCli(raw);
    std::vector<std::string> argv{"riptide_sim_predict"};
    argv.insert(argv.end(), args.begin(), args.end());
    const auto result = cli.parse(lyra::args(argv.begin(), argv.end()));
    if (!result) {
        throw CliError(result.message());
    }

    PredictCliOptions opts;
    opts.show_help = raw.show_help;
    opts.show_version = raw.show_version;
    if (opts.show_help || opts.show_version) {
        return opts;
    }
    if (raw.config_path.empty()) {
        throw CliError("--config is required");
    }
    opts.config_path = raw.config_path;
    if (raw.events.empty()) {
        throw CliError("--events is required");
    }
    opts.events = parseCliNumber<std::uint64_t>(raw.events, "--events");
    if (opts.events == 0) {
        throw CliError("--events: must be positive");
    }
    opts.bins_u = parseBins(raw.bins_u, "--bins-u");
    opts.bins_z = parseBins(raw.bins_z, "--bins-z");
    if (!raw.energy_nodes.empty()) {
        opts.energy_nodes = parseCliNumber<std::size_t>(raw.energy_nodes, "--nodes");
        if (opts.energy_nodes == 0 || opts.energy_nodes > kMaxEnergyNodes) {
            throw CliError("--nodes: must be in [1, " + std::to_string(kMaxEnergyNodes) + "]");
        }
    }
    if (!raw.data_dir.empty()) {
        opts.data_dir = raw.data_dir;
    }
    return opts;
}

std::string predictHelpText() {
    RawPredictOptions raw;
    std::ostringstream out;
    out << makePredictCli(raw);
    return out.str();
}

std::string checkHelpText() {
    RawOptions raw;
    std::ostringstream out;
    out << makeCli(raw);
    return out.str();
}

} // namespace riptide::sim::checks
