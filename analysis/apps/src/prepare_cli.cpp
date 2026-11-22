#include "riptide/apps/prepare_cli.hpp"

#include <lyra/lyra.hpp>

#include <cmath>
#include <exception>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "riptide/common/data_paths.hpp"

namespace riptide::apps {

namespace {

struct RawOptions {
    bool show_help = false;
    bool show_version = false;
    std::string data_dir;
    std::string run_id;
    std::string analysis_id;
    std::string seed;
    std::string sigma_Ep;
    std::string sigma_theta;
    std::string cuts;
    std::string hp_check_id;
    std::string max_events;
};

lyra::cli makeCli(RawOptions& raw) {
    return lyra::cli() | lyra::help(raw.show_help) |
           lyra::opt(raw.show_version)["--version"]("print the versions and exit") |
           lyra::opt(raw.run_id, "id")["--run-id"]("run to read, under <root>/sim/, required") |
           lyra::opt(raw.analysis_id, "id")["--analysis-id"](
               "output directory name under <root>/analysis/<run-id>/, required; must not "
               "exist") |
           lyra::opt(raw.seed, "n")["--seed"]("seed of the resolution draws, required") |
           lyra::opt(raw.sigma_Ep,
                     "MeV")["--sigma-Ep"]("Gaussian resolution of Ep in MeV, >= 0, required") |
           lyra::opt(raw.sigma_theta, "rad")["--sigma-theta"](
               "Gaussian resolution of theta_p and of each tangent component of the track "
               "direction in rad, >= 0, required") |
           lyra::opt(raw.cuts, "list")["--cuts"](
               "comma-separated cut chain, e.g. single_np,contained,evis_min=0.5 (MeV), "
               "required") |
           lyra::opt(raw.hp_check_id, "id")["--hp-from-check"](
               "take the HP target mass from <root>/validation/<run-id>/<id>/checks.json; "
               "required for run.json version 1") |
           lyra::opt(raw.max_events, "n")["--max-events"](
               "read only the first n events by event_id (default: all)") |
           lyra::opt(raw.data_dir, "dir")["--data-dir"]("data root (default: $RIPTIDE_DATA_DIR, "
                                                        "then /mnt/external_ssd/riptide-sim)");
}

std::string requireValue(const std::string& value, std::string_view option) {
    if (value.empty()) {
        throw CliError(std::string{option} + " is required");
    }
    return value;
}

std::string requireId(const std::string& value, std::string_view option) {
    requireValue(value, option);
    if (!data::isValidId(value)) {
        throw CliError(std::string{option} +
                       ": only [A-Za-z0-9._-], not starting with '.' or '-', at most " +
                       std::to_string(data::kMaxIdLength) + " characters");
    }
    return value;
}

double parseSigma(const std::string& value, std::string_view option) {
    const auto sigma = parseCliNumber<double>(requireValue(value, option), option);
    // from_chars accepts "inf" and "nan"; neither is a resolution.
    if (!std::isfinite(sigma) || sigma < 0.0) {
        throw CliError(std::string{option} + ": must be finite and >= 0");
    }
    return sigma;
}

} // namespace

PrepareOptions parsePrepareCommandLine(const std::vector<std::string>& args) {
    RawOptions raw;
    const auto cli = makeCli(raw);
    std::vector<std::string> argv{"riptide_analysis_prepare"};
    argv.insert(argv.end(), args.begin(), args.end());
    const auto result = cli.parse(lyra::args(argv.begin(), argv.end()));
    if (!result) {
        throw CliError(result.message());
    }

    PrepareOptions opts;
    opts.show_help = raw.show_help;
    opts.show_version = raw.show_version;
    if (opts.show_help || opts.show_version) {
        return opts;
    }
    opts.run_id = requireId(raw.run_id, "--run-id");
    opts.analysis_id = requireId(raw.analysis_id, "--analysis-id");
    opts.seed = parseCliNumber<std::uint64_t>(requireValue(raw.seed, "--seed"), "--seed");
    opts.resolution.sigma_Ep = parseSigma(raw.sigma_Ep, "--sigma-Ep");
    opts.resolution.sigma_theta = parseSigma(raw.sigma_theta, "--sigma-theta");
    try {
        opts.selection = io::Selection::parse(requireValue(raw.cuts, "--cuts"));
    } catch (const CliError&) {
        throw;
    } catch (const std::exception& e) {
        throw CliError(std::string{"--cuts: "} + e.what());
    }
    if (!raw.hp_check_id.empty()) {
        opts.hp_check_id = requireId(raw.hp_check_id, "--hp-from-check");
    }
    if (!raw.max_events.empty()) {
        opts.max_events = parseCliNumber<std::uint64_t>(raw.max_events, "--max-events");
        if (opts.max_events == 0) {
            throw CliError("--max-events: must be positive");
        }
    }
    if (!raw.data_dir.empty()) {
        opts.data_dir = raw.data_dir;
    }
    return opts;
}

std::string prepareHelpText() {
    RawOptions raw;
    std::ostringstream out;
    out << makeCli(raw);
    return out.str();
}

} // namespace riptide::apps
