#pragma once

// Run configuration of the simulation, read from a versioned JSON file in
// config/. Parsing and validation are pure functions with no Geant4
// dependency, so they are unit-tested in the Geant4-free CI build.
//
// JSON keys carry their unit as a suffix (_mm, _MeV). Values are converted to
// the project units (riptide::units) on the way in, so the struct below holds
// plain doubles in MeV and mm.

#include <nlohmann/json_fwd.hpp>

#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>

#include "riptide/common/units.hpp"

namespace riptide::sim {

// Version of the JSON layout below. New runs need exactly this version
// (loadRunConfig); parseRunConfig also reads the older layouts still found in
// the run.json of existing runs, from kOldestConfigVersion on.
inline constexpr int kConfigVersion = 3;
inline constexpr int kOldestConfigVersion = 1;

// Version 2 added the Birks constant and the step function of the scintillator
// physics; a version 1 file has neither.
inline constexpr int kFirstBirksConfigVersion = 2;

// Version 3 added source.beam (pencil or far-field); an older file always
// describes a pencil beam.
inline constexpr int kFirstBeamConfigVersion = 3;

// Unit of the Birks constant in the JSON file (key suffix _mm_per_MeV).
inline constexpr double kBirksUnit = units::mm / units::MeV;

// Bounds accepted by Geant4 (G4EmExtraParameters::SetStepFunctionMuHad), which
// otherwise ignores the values with only a warning: checking them here turns
// that warning into an error before the run starts.
inline constexpr double kMaxDRoverRange = 1.0;

// G4NeutronHP (used by the *_HP physics lists) relies on evaluated neutron data
// that end at 20 MeV; above that the model changes and the n-p truth would no
// longer come from the evaluated cross sections.
inline constexpr double kNeutronHpMaxEnergy = 20.0 * units::MeV;

// Tolerance on |direction| - 1 in the input file: the direction is normalised
// after the check, so this only catches typos, not rounding.
inline constexpr double kDirectionNormTolerance = 1.0e-6;

class ConfigError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

enum class SourceMode { Mono, Flat };

// Pencil: every neutron starts at one point along one direction. FarField: a
// parallel beam from a distant source, uniform over a disk perpendicular to the
// direction that covers the whole scintillator (source_sampling.hpp).
enum class BeamMode { Pencil, FarField };

struct GeometryConfig {
    Vec3 scint_size; // full edge lengths of the scintillator box, centred at the origin [mm]
};

struct SourceConfig {
    BeamMode beam = BeamMode::Pencil;
    SourceMode mode = SourceMode::Mono;
    double En_min = 0.0;    // [MeV]; equal to En_max for a monoenergetic beam
    double En_max = 0.0;    // [MeV]
    bool isotropic = false; // far-field only: a new direction per event, uniform on the sphere
    Vec3 direction;         // unit vector; unused when isotropic
    Vec3 position;          // pencil only: gun position [mm], outside the scintillator
    double margin = 0.0;    // far-field only: start-disk distance beyond the half-diagonal [mm]
};

struct ScintillatorConfig {
    double birks_kB = 0.0;    // Birks constant of the scintillator [mm/MeV]; 0 disables quenching
    std::string birks_source; // where the value comes from (a decision id), recorded in run.json
};

// Geant4 step function of the ionisation of muons and hadrons, protons
// included: a step is at most max(dRoverRange * range, finalRange). Birks is
// applied step by step, so the visible energy of a stopping proton depends on
// it and the values are part of the run configuration.
struct StepFunction {
    double dRoverRange = 0.0; // dimensionless, in (0, 1]
    double finalRange = 0.0;  // [mm], positive
};

struct RunConfig {
    int config_version = kConfigVersion;
    GeometryConfig geometry;
    SourceConfig source;
    std::string physics_list;
    // Present from kFirstBirksConfigVersion on; empty only for a version 1
    // configuration read back from the run.json of an older run.
    std::optional<ScintillatorConfig> scintillator;
    std::optional<StepFunction> proton_step_function;
};

// Half edge lengths of the scintillator box [mm].
[[nodiscard]] Vec3 scintHalfSize(const GeometryConfig& geometry);

// Half-size of the cubic vacuum world: it must contain both the scintillator
// and every start point of the source, with room to spare so no primary sits
// on the world boundary (a primary starting on the boundary is killed at the
// first step).
inline constexpr double kWorldScale = 1.5;
inline constexpr double kWorldMargin = 10.0 * units::mm;
[[nodiscard]] double worldHalfSize(const RunConfig& config);

// Throws ConfigError with the offending key on any missing, unknown or invalid
// field: a typo in a key must never silently fall back to a default. Accepts
// every version from kOldestConfigVersion to kConfigVersion.
[[nodiscard]] RunConfig parseRunConfig(const nlohmann::json& doc);

// Reads and parses the configuration of a new run; throws ConfigError on I/O
// or syntax errors and on any version other than kConfigVersion, because an
// older file lacks physics settings that a new run must state explicitly.
[[nodiscard]] RunConfig loadRunConfig(const std::filesystem::path& path);

// Inverse of parseRunConfig, used to record the effective configuration in
// run.json: parseRunConfig(toJson(c)) reproduces c.
[[nodiscard]] nlohmann::json toJson(const RunConfig& config);

} // namespace riptide::sim
