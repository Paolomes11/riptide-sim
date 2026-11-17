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
//
// Both seeds are also written to every event as signed 32-bit columns (the
// Geant4 ntuple has no 64-bit integer column), so they are capped at 2^31 - 1.

#include <cstdint>

namespace riptide::sim {

// MixMaxRng::setSeeds ignores the bits above 32, and the seed columns are
// signed 32-bit: the run seed must fit in both.
inline constexpr std::uint64_t kMaxSeed = 0x7FFFFFFFULL;
// Event ids are G4int (BeamOn takes a G4int), so they always fit in 32 bits.
inline constexpr std::int64_t kMaxSeededEvents = 0x7FFFFFFFLL;

struct EventSeeds {
    std::int32_t seed_0 = 0; // first seed passed to the engine: the event id
    std::int32_t seed_1 = 0; // second seed: the run seed
};

// Preconditions (checked by the callers): 0 <= event_id <= kMaxSeededEvents and
// run_seed <= kMaxSeed; the values then fit unchanged in 31 bits.
[[nodiscard]] constexpr EventSeeds eventSeeds(std::uint64_t run_seed, std::int64_t event_id) {
    return {static_cast<std::int32_t>(event_id), static_cast<std::int32_t>(run_seed)};
}

} // namespace riptide::sim
