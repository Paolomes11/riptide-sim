#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <set>
#include <utility>

#include "riptide/sim/cli_options.hpp"
#include "riptide/sim/event_seeds.hpp"

namespace {

using riptide::sim::eventSeeds;
using riptide::sim::kMaxSeed;
using riptide::sim::kMaxSeededEvents;

constexpr std::uint64_t kMask32 = 0xFFFFFFFFULL;

// What MixMaxRng::setSeeds actually uses of a seed.
std::uint64_t low32(std::int64_t seed) {
    return static_cast<std::uint64_t>(seed) & kMask32;
}

} // namespace

TEST_CASE("event seeds are the event id and the run seed", "[seeds]") {
    const auto s = eventSeeds(riptide::sim::kDefaultSeed, 42);
    CHECK(s.seed_0 == 42);
    CHECK(s.seed_1 == 20260907);
}

TEST_CASE("event seeds survive the engine's 32-bit truncation", "[seeds]") {
    // At the limits the values the engine sees equal the values recorded, so
    // the recorded pair reproduces the event and no two pairs alias.
    const auto s = eventSeeds(kMaxSeed, kMaxSeededEvents);
    CHECK(low32(s.seed_0) == static_cast<std::uint64_t>(kMaxSeededEvents));
    CHECK(low32(s.seed_1) == kMaxSeed);
    CHECK(s.seed_0 >= 0);
    CHECK(s.seed_1 >= 0);
}

TEST_CASE("distinct events or runs give distinct engine seeds", "[seeds]") {
    std::set<std::pair<std::uint64_t, std::uint64_t>> seen;
    for (const std::uint64_t run_seed : {std::uint64_t{0}, std::uint64_t{1}, kMaxSeed}) {
        for (std::int64_t id = 0; id < 100; ++id) {
            const auto s = eventSeeds(run_seed, id);
            CHECK(seen.emplace(low32(s.seed_0), low32(s.seed_1)).second);
        }
    }
}

TEST_CASE("the CLI limits keep seeds within 32 bits", "[seeds]") {
    STATIC_REQUIRE(riptide::sim::kMaxEvents - 1 <= static_cast<std::uint64_t>(kMaxSeededEvents));
    STATIC_REQUIRE(riptide::sim::kDefaultSeed <= kMaxSeed);
}
