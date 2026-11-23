#pragma once

// The pure part of riptide_case_a_geant4 (D37, D38): turning selected Geant4
// events into the sample the Case A engine reads, the --truth spec of the toy
// source, and how many disjoint experiments a row can have. None of it reads
// files, so it is built and tested in every preset.

#include <cstddef>
#include <string_view>

#include "riptide/apps/prepare_cli.hpp"
#include "riptide/common/event_sample.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/stats/pseudo_experiment.hpp"

namespace riptide::apps {

// The pilot divides every count by this, keeping at least kPilotMinEntries so
// an ensemble keeps its error; N per experiment is never scaled.
inline constexpr std::size_t kPilotDivisor = 10;
inline constexpr std::size_t kPilotMinEntries = 2;

// An ensemble summary needs at least two entries for its errors.
inline constexpr std::size_t kMinExperiments = 2;

// Relative tolerance on truth_En of a mono run: every event carries the
// configured energy, copied as a double, so only rounding may differ.
inline constexpr double kMonoTruthRelTolerance = 1e-12;

[[nodiscard]] std::size_t pilotScaled(std::size_t count, bool pilot) noexcept;

// The columns the Case A engine reads. truth_Ep is the physical proton energy
// (truth_Ep_phys): the engine's likelihood is written for n-p kinematics
// without the HP target-mass convention. Throws std::invalid_argument unless
// truth_En, truth_Ep_phys, truth_theta_p, reco_Ep and reco_theta_p all have
// one entry per event.
[[nodiscard]] stats::PseudoSample pseudoSampleFromEvents(const EventSample& events);

// Parses "mono:<En>" or "flat:<min>:<max>" (MeV). Throws CliError on anything
// else, on a value that is not positive and finite, or on min >= max.
[[nodiscard]] io::SourceSpectrum parseTruth(std::string_view spec);

// Experiments of a row: the planned M (pilot-scaled) or fewer when the pool
// has not enough events for M disjoint experiments of N events. Returns 0 if
// fewer than kMinExperiments fit, and the row is then skipped.
[[nodiscard]] std::size_t rowExperiments(std::size_t planned, std::size_t available_events,
                                         std::size_t events_per_experiment, bool pilot);

// Throws std::invalid_argument unless every truth_En equals En within
// kMonoTruthRelTolerance: a combined posterior assumes one common E_n.
void requireMonoTruth(const stats::PseudoSample& sample, double En);

} // namespace riptide::apps
