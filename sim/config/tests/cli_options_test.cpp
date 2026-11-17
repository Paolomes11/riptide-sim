#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <string>
#include <vector>

#include "riptide/sim/cli_options.hpp"
#include "riptide/sim/event_seeds.hpp"

namespace {

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using riptide::sim::CliError;
using riptide::sim::parseCommandLine;

using Args = std::vector<std::string>;

Args required() {
    return {"--config", "config/mono_14MeV.json", "--events", "1000", "--run-id", "smoke_001"};
}

Args with(Args args, const Args& extra) {
    args.insert(args.end(), extra.begin(), extra.end());
    return args;
}

void requireError(const Args& args, const std::string& fragment) {
    REQUIRE_THROWS_MATCHES(parseCommandLine(args), CliError,
                           MessageMatches(ContainsSubstring(fragment)));
}

} // namespace

TEST_CASE("required options alone give the documented defaults", "[cli]") {
    const auto o = parseCommandLine(required());
    CHECK_FALSE(o.show_help);
    CHECK_FALSE(o.show_version);
    CHECK(o.config == "config/mono_14MeV.json");
    CHECK(o.events == 1000);
    CHECK(o.run_id == "smoke_001");
    CHECK(o.threads == 1);
    CHECK(o.seed == riptide::sim::kDefaultSeed);
    CHECK(o.seed == 20260907U);
    CHECK_FALSE(o.data_dir.has_value());
    CHECK_FALSE(o.allow_system_disk);
}

TEST_CASE("optional options are read", "[cli]") {
    const auto o =
        parseCommandLine(with(required(), {"--threads", "4", "--seed", "42", "--data-dir",
                                           "/data/x", "--allow-system-disk"}));
    CHECK(o.threads == 4);
    CHECK(o.seed == 42U);
    REQUIRE(o.data_dir.has_value());
    CHECK(*o.data_dir == "/data/x");
    CHECK(o.allow_system_disk);
}

TEST_CASE("the seed is limited to 32 bits", "[cli]") {
    const auto o = parseCommandLine(with(required(), {"--seed", "4294967295"}));
    CHECK(o.seed == riptide::sim::kMaxSeed);
    requireError(with(required(), {"--seed", "4294967296"}), "--seed");
    requireError(with(required(), {"--seed", "18446744073709551616"}), "--seed");
}

TEST_CASE("--help and --version need no other option", "[cli]") {
    CHECK(parseCommandLine({"--help"}).show_help);
    CHECK(parseCommandLine({"--version"}).show_version);
    CHECK_THAT(riptide::sim::helpText(), ContainsSubstring("--config"));
}

TEST_CASE("missing required options are reported", "[cli]") {
    requireError({"--events", "1", "--run-id", "a"}, "--config");
    requireError({"--config", "c.json", "--run-id", "a"}, "--events");
    requireError({"--config", "c.json", "--events", "1"}, "--run-id");
}

TEST_CASE("malformed values are rejected", "[cli]") {
    requireError(with(required(), {"--threads", "0"}), "--threads");
    requireError(with(required(), {"--threads", "100000"}), "--threads");
    requireError(with(required(), {"--threads", "2x"}), "--threads");
    requireError(with(required(), {"--seed", "abc"}), "--seed");
    requireError({"--config", "c.json", "--events", "0", "--run-id", "a"}, "--events");
    requireError({"--config", "c.json", "--events", "1e3", "--run-id", "a"}, "--events");
    requireError({"--config", "c.json", "--events", "2147483648", "--run-id", "a"}, "--events");
    requireError({"--config", "c.json", "--events", "1", "--run-id", "../evil"}, "--run-id");
}

TEST_CASE("negative numbers never wrap to large unsigned values", "[cli]") {
    // Either Lyra refuses "-5" as a value or the conversion does; both throw.
    CHECK_THROWS_AS(parseCommandLine(with(required(), {"--seed", "-5"})), CliError);
    CHECK_THROWS_AS(parseCommandLine({"--config", "c.json", "--events", "-1", "--run-id", "a"}),
                    CliError);
}

TEST_CASE("unknown options are rejected", "[cli]") {
    CHECK_THROWS_AS(parseCommandLine(with(required(), {"--evnets", "10"})), CliError);
}
