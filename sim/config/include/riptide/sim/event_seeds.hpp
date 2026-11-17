#pragma once

// Per-event seeding of the Geant4 random engine.
//
// Every event is seeded from (event_id, run seed) alone, at the start of its
// primary generation, so an event is reproducible on its own and the output does
// not depend on the number of threads or on which thread ran the event.
//
// The pair is handed to CLHEP::MixMaxRng::setSeeds, which keeps the low 32 bits
// of each seed and derives the event stream from the IDs with MixMax's skip-ahead
// (seed_uniquestream). The MixMax authors guarantee that streams derived from IDs
// differing in at least one bit do not overlap (within 1e100 draws), so distinct
// (event_id, run seed) pairs give non-overlapping streams by construction, with
// no hashing and no collision probability. This holds only if both values fit in
// 32 bits, hence kMaxSeed and kMaxSeededEvents.

#include <cstdint>

namespace riptide::sim {

// MixMaxRng::setSeeds ignores the bits above 32: larger seeds would alias.
inline constexpr std::uint64_t kMaxSeed = 0xFFFFFFFFULL;
// Event ids are G4int (BeamOn takes a G4int), so they always fit in 32 bits.
inline constexpr std::int64_t kMaxSeededEvents = 0x7FFFFFFFLL;

struct EventSeeds {
    std::int64_t seed_0 = 0; // first seed passed to the engine: the event id
    std::int64_t seed_1 = 0; // second seed: the run seed
};

// Preconditions (checked by the callers): 0 <= event_id <= kMaxSeededEvents and
// run_seed <= kMaxSeed; the values then fit unchanged in the engine's 32 bits.
[[nodiscard]] constexpr EventSeeds eventSeeds(std::uint64_t run_seed, std::int64_t event_id) {
    return {event_id, static_cast<std::int64_t>(run_seed)};
}

} // namespace riptide::sim
