#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <string>
#include <vector>

#include "riptide/sim/check_cli_options.hpp"

namespace {

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using riptide::sim::CliError;
using riptide::sim::checks::kDefaultBins;
using riptide::sim::checks::kMaxBins;
using riptide::sim::checks::parseCheckCommandLine;
using riptide::sim::checks::parsePredictCommandLine;

using Args = std::vector<std::string>;

Args required() {
    return {"--run-id", "mono_14MeV_1e7", "--check-id", "first"};
}

Args with(Args args, const Args& extra) {
    args.insert(args.end(), extra.begin(), extra.end());
    return args;
}

void requireError(const Args& args, const std::string& fragment) {
    REQUIRE_THROWS_MATCHES(parseCheckCommandLine(args), CliError,
                           MessageMatches(ContainsSubstring(fragment)));
}

void requirePredictError(const Args& args, const std::string& fragment) {
    REQUIRE_THROWS_MATCHES(parsePredictCommandLine(args), CliError,
                           MessageMatches(ContainsSubstring(fragment)));
}

} // namespace

TEST_CASE("check cli: required options alone give the defaults", "[check_cli]") {
    const auto o = parseCheckCommandLine(required());
    CHECK_FALSE(o.show_help);
    CHECK_FALSE(o.show_version);
    CHECK(o.run_id == "mono_14MeV_1e7");
    CHECK(o.check_id == "first");
    CHECK(o.bins_u == kDefaultBins);
    CHECK(o.bins_z == kDefaultBins);
    CHECK_FALSE(o.data_dir.has_value());
}

TEST_CASE("check cli: optional options are read", "[check_cli]") {
    const auto o = parseCheckCommandLine(
        with(required(), {"--bins-u", "50", "--bins-z", "1", "--data-dir", "/data/x"}));
    CHECK(o.bins_u == 50);
    CHECK(o.bins_z == 1);
    REQUIRE(o.data_dir.has_value());
    CHECK(*o.data_dir == "/data/x");
}

TEST_CASE("check cli: bin counts are bounded", "[check_cli]") {
    CHECK(parseCheckCommandLine(with(required(), {"--bins-u", std::to_string(kMaxBins)})).bins_u ==
          kMaxBins);
    requireError(with(required(), {"--bins-u", "0"}), "--bins-u");
    requireError(with(required(), {"--bins-z", std::to_string(kMaxBins + 1)}), "--bins-z");
    requireError(with(required(), {"--bins-u", "-3"}), "--bins-u");
    requireError(with(required(), {"--bins-z", "20x"}), "--bins-z");
}

TEST_CASE("check cli: ids are required and validated", "[check_cli]") {
    requireError({"--check-id", "a"}, "--run-id");
    requireError({"--run-id", "a"}, "--check-id");
    requireError({"--run-id", "../x", "--check-id", "a"}, "--run-id");
    requireError({"--run-id", "a", "--check-id", ".hidden"}, "--check-id");
}

TEST_CASE("check cli: --help and --version need no other option", "[check_cli]") {
    CHECK(parseCheckCommandLine({"--help"}).show_help);
    CHECK(parseCheckCommandLine({"--version"}).show_version);
    const std::string help = riptide::sim::checks::checkHelpText();
    CHECK_THAT(help, ContainsSubstring("--check-id"));
    CHECK_THAT(help, ContainsSubstring("default " + std::to_string(kDefaultBins)));
}

TEST_CASE("check cli: unknown options are rejected", "[check_cli]") {
    REQUIRE_THROWS_AS(parseCheckCommandLine(with(required(), {"--events", "3"})), CliError);
}

TEST_CASE("predict cli: required options alone give the defaults", "[check_cli]") {
    const auto o =
        parsePredictCommandLine({"--config", "config/mono_14MeV.json", "--events", "10000000"});
    CHECK(o.config_path == "config/mono_14MeV.json");
    CHECK(o.events == 10000000);
    CHECK(o.bins_u == kDefaultBins);
    CHECK(o.bins_z == kDefaultBins);
    CHECK(o.energy_nodes == riptide::sim::checks::kDefaultEnergyNodes);
    CHECK_FALSE(o.data_dir.has_value());
}

TEST_CASE("predict cli: optional options are read and bounded", "[check_cli]") {
    const Args base{"--config", "c.json", "--events", "5"};
    const auto o = parsePredictCommandLine(
        with(base, {"--bins-u", "50", "--bins-z", "40", "--nodes", "128", "--data-dir", "/d"}));
    CHECK(o.bins_u == 50);
    CHECK(o.bins_z == 40);
    CHECK(o.energy_nodes == 128);
    REQUIRE(o.data_dir.has_value());
    CHECK(*o.data_dir == "/d");
    requirePredictError(with(base, {"--nodes", "0"}), "--nodes");
    requirePredictError(
        with(base, {"--nodes", std::to_string(riptide::sim::checks::kMaxEnergyNodes + 1)}),
        "--nodes");
    requirePredictError(with(base, {"--bins-u", "0"}), "--bins-u");
}

TEST_CASE("predict cli: config and a positive event count are required", "[check_cli]") {
    requirePredictError({"--events", "5"}, "--config");
    requirePredictError({"--config", "c.json"}, "--events");
    requirePredictError({"--config", "c.json", "--events", "0"}, "--events");
    requirePredictError({"--config", "c.json", "--events", "1e7"}, "--events");
    CHECK(parsePredictCommandLine({"--help"}).show_help);
    CHECK(parsePredictCommandLine({"--version"}).show_version);
    CHECK_THAT(riptide::sim::checks::predictHelpText(), ContainsSubstring("--nodes"));
    REQUIRE_THROWS_AS(parsePredictCommandLine({"--config", "c", "--events", "5", "--run-id", "x"}),
                      CliError);
}
