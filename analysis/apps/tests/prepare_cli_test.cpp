#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <string>
#include <vector>

#include "riptide/apps/prepare_cli.hpp"
#include "riptide/io/selection.hpp"

namespace {

using riptide::apps::CliError;
using riptide::apps::parsePrepareCommandLine;
using Args = std::vector<std::string>;

Args required() {
    return {"--run-id",   "r1",  "--analysis-id", "a1",   "--seed", "20260907",
            "--sigma-Ep", "0.1", "--sigma-theta", "0.08", "--cuts", "single_np,contained"};
}

Args with(Args args, const Args& extra) {
    args.insert(args.end(), extra.begin(), extra.end());
    return args;
}

// The required arguments without the option `name` and its value.
Args without(const std::string& name) {
    Args args = required();
    for (auto it = args.begin(); it != args.end(); ++it) {
        if (*it == name) {
            args.erase(it, it + 2);
            break;
        }
    }
    return args;
}

} // namespace

TEST_CASE("the required options are parsed", "[prepare_cli]") {
    const auto opts = parsePrepareCommandLine(required());
    CHECK(opts.run_id == "r1");
    CHECK(opts.analysis_id == "a1");
    CHECK(opts.seed == 20260907U);
    CHECK_THAT(opts.resolution.sigma_Ep, Catch::Matchers::WithinULP(0.1, 0));
    CHECK_THAT(opts.resolution.sigma_theta, Catch::Matchers::WithinULP(0.08, 0));
    CHECK(opts.selection.toString() == "single_np,contained");
    CHECK_FALSE(opts.hp_check_id.has_value());
    CHECK(opts.max_events == 0U);
    CHECK_FALSE(opts.data_dir.has_value());
}

TEST_CASE("the optional options are parsed", "[prepare_cli]") {
    const auto opts = parsePrepareCommandLine(with(
        required(), {"--hp-from-check", "c1", "--max-events", "100000", "--data-dir", "/x/y"}));
    CHECK(opts.hp_check_id == "c1");
    CHECK(opts.max_events == 100000U);
    CHECK(opts.data_dir == "/x/y");
}

TEST_CASE("zero resolutions are accepted", "[prepare_cli]") {
    const auto opts =
        parsePrepareCommandLine({"--run-id", "r1", "--analysis-id", "a1", "--seed", "1",
                                 "--sigma-Ep", "0", "--sigma-theta", "0", "--cuts", "single_np"});
    CHECK_THAT(opts.resolution.sigma_Ep, Catch::Matchers::WithinAbs(0.0, 0.0));
    CHECK_THAT(opts.resolution.sigma_theta, Catch::Matchers::WithinAbs(0.0, 0.0));
}

TEST_CASE("every analysis setting is required", "[prepare_cli]") {
    for (const char* name :
         {"--run-id", "--analysis-id", "--seed", "--sigma-Ep", "--sigma-theta", "--cuts"}) {
        INFO(name);
        CHECK_THROWS_AS(parsePrepareCommandLine(without(name)), CliError);
    }
}

TEST_CASE("malformed values are refused", "[prepare_cli]") {
    const auto replaced = [](const std::string& name, const std::string& value) {
        return with(without(name), {name, value});
    };
    for (const auto& [name, value] : std::vector<std::pair<std::string, std::string>>{
             {"--seed", "-1"},
             {"--seed", "12x"},
             {"--seed", "1e3"},
             {"--sigma-Ep", "-0.1"},
             {"--sigma-Ep", "nan"},
             {"--sigma-Ep", "inf"},
             {"--sigma-theta", "0.08rad"},
             {"--sigma-theta", ""},
             {"--cuts", "single_np,unknown"},
             {"--cuts", "evis_min"},
             {"--cuts", "single_np,single_np"},
             {"--run-id", "../r1"},
             {"--analysis-id", "-a"},
         }) {
        INFO(name << " = '" << value << "'");
        CHECK_THROWS_AS(parsePrepareCommandLine(replaced(name, value)), CliError);
    }
    CHECK_THROWS_AS(parsePrepareCommandLine(with(required(), {"--max-events", "0"})), CliError);
    CHECK_THROWS_AS(parsePrepareCommandLine(with(required(), {"--hp-from-check", "a/b"})),
                    CliError);
    CHECK_THROWS_AS(parsePrepareCommandLine(with(required(), {"--unknown"})), CliError);
}

TEST_CASE("help and version need no other option", "[prepare_cli]") {
    CHECK(parsePrepareCommandLine({"--help"}).show_help);
    CHECK(parsePrepareCommandLine({"--version"}).show_version);
    const auto help = riptide::apps::prepareHelpText();
    for (const char* name : {"--run-id", "--analysis-id", "--seed", "--sigma-Ep", "--sigma-theta",
                             "--cuts", "--hp-from-check", "--max-events", "--data-dir"}) {
        INFO(name);
        CHECK(help.find(name) != std::string::npos);
    }
}
