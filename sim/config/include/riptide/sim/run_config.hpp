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
#include <stdexcept>
#include <string>

#include "riptide/common/units.hpp"

namespace riptide::sim {

// Version of the JSON layout below; files with another version are refused.
inline constexpr int kConfigVersion = 1;

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

struct GeometryConfig {
    Vec3 scint_size; // full edge lengths of the scintillator box, centred at the origin [mm]
};

struct SourceConfig {
    SourceMode mode = SourceMode::Mono;
    double En_min = 0.0; // [MeV]; equal to En_max for a monoenergetic beam
    double En_max = 0.0; // [MeV]
    Vec3 direction;      // unit vector
    Vec3 position;       // gun position [mm], outside the scintillator
};

struct RunConfig {
    int config_version = kConfigVersion;
    GeometryConfig geometry;
    SourceConfig source;
    std::string physics_list;
};

// Throws ConfigError with the offending key on any missing, unknown or invalid
// field: a typo in a key must never silently fall back to a default.
[[nodiscard]] RunConfig parseRunConfig(const nlohmann::json& doc);

// Reads and parses a JSON file; throws ConfigError on I/O or syntax errors.
[[nodiscard]] RunConfig loadRunConfig(const std::filesystem::path& path);

// Inverse of parseRunConfig, used to record the effective configuration in
// run.json: parseRunConfig(toJson(c)) reproduces c.
[[nodiscard]] nlohmann::json toJson(const RunConfig& config);

} // namespace riptide::sim
