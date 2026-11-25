#pragma once

// The pure part of riptide_case_b_geant4 (D45): turning selected Geant4 events
// into the sample the Case B engine reads, the fixed neutron direction of a
// run and of the --direction option, and the check that every event carries
// it. None of it reads files, so it is built and tested in every preset.

#include <nlohmann/json.hpp>

#include <string_view>

#include "riptide/apps/prepare_cli.hpp"
#include "riptide/common/direction.hpp"
#include "riptide/common/event_sample.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"

namespace riptide::apps {

// Tolerance on |direction| - 1 of a configured direction, the one of the
// simulation config: the direction is normalised after the check, so this
// only catches typos, not rounding.
inline constexpr double kDirectionNormTolerance = 1e-6;

// Tolerance on |track| - 1 of a proton direction in the events: the
// simulation writes unit momentum directions and the response keeps the
// smeared one unit, so only rounding may differ.
inline constexpr double kTrackNormTolerance = 1e-9;

// Angle between truth_omega_n of an event and the configured direction: the
// simulation copies the normalised direction into every event, so only
// rounding may differ.
inline constexpr double kFixedDirectionTolerance = 1e-9; // rad

// The Case B sample of selected events: truth_Ep_phys as the truth Ep, the
// truth and reco proton directions as tracks. Throws std::invalid_argument
// if a column does not have one entry per event or a track is not a unit
// vector within kTrackNormTolerance (an event without a proton has none).
[[nodiscard]] stats::CaseBPseudoSample caseBSampleFromEvents(const EventSample& events);

// "x,y,z" as a unit vector, normalised after a norm check within
// kDirectionNormTolerance. Throws CliError on a malformed spec.
[[nodiscard]] geometry::Direction parseDirection(std::string_view spec);

// The fixed neutron direction of a run.json document, config.source.direction,
// normalised as the simulation does. Throws io::IoError if it is missing,
// "isotropic", not three finite numbers or not a unit vector within
// kDirectionNormTolerance.
[[nodiscard]] geometry::Direction fixedDirectionFromRunRecord(const nlohmann::json& record);

// Throws std::invalid_argument unless every event has truth_omega_n within
// kFixedDirectionTolerance of omega_n: a Case B experiment assumes one
// direction shared by all its events.
void requireFixedDirection(const EventSample& events, const geometry::Direction& omega_n);

} // namespace riptide::apps
