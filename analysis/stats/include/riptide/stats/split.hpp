#pragma once

// Disjoint experiments drawn from one finite sample of selected events, as
// for a Geant4 run: no event enters two experiments of the same row, so the
// experiments are independent and the binomial error of a coverage is
// honest. The split depends only on the event count and the engine, never
// on the event values.

#include <cstddef>
#include <vector>

#include "riptide/common/random.hpp"
#include "riptide/stats/pseudo_experiment.hpp"

namespace riptide::stats {

// `experiments` blocks of `events_per_experiment` indices each, taken from
// 0 .. events - 1. A partial Fisher-Yates shuffle orders the first
// experiments * events_per_experiment positions, which is a uniform random
// ordered subset; the blocks are consecutive slices of it and the leftover
// events are unused. Throws if events_per_experiment == 0 or the blocks do
// not fit; zero experiments give an empty result.
[[nodiscard]] std::vector<std::vector<std::size_t>>
disjointExperiments(std::size_t events, std::size_t events_per_experiment, std::size_t experiments,
                    random::Engine& engine);

// Largest number of disjoint experiments of events_per_experiment events
// that fit in `events`. Throws if events_per_experiment == 0.
[[nodiscard]] std::size_t maxDisjointExperiments(std::size_t events,
                                                 std::size_t events_per_experiment);

// The events at `indices`, in that order. Throws if the columns of `sample`
// have different lengths or an index is out of range.
[[nodiscard]] PseudoSample pseudoSampleSubset(const PseudoSample& sample,
                                              const std::vector<std::size_t>& indices);

} // namespace riptide::stats
