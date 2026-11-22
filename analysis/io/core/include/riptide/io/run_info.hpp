#pragma once

// What the analysis needs to know about a simulation run, read from its
// run.json and, for the HP target mass of runs written before run.json
// carried it, from checks.json of a check of the same run.
//
// The analysis does not link the simulation: the layout of both files is the
// contract written in docs/architecture/data_layout.md, and this reader
// refuses what it does not know instead of guessing.

#include <nlohmann/json.hpp>

#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace riptide::io {

class IoError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// run.json layouts this reader understands. Version 2 added hp_elastic_H1.
inline constexpr int kOldestRunRecordVersion = 1;
inline constexpr int kNewestRunRecordVersion = 2;
inline constexpr int kFirstRunRecordWithHp = 2;

// target_mass_MeV is written as awr * kNeutronMass in double precision: the two
// agree to the last bits unless one of them was edited or mis-copied.
inline constexpr double kHpMassRelTolerance = 1e-12;

inline constexpr std::string_view kRunRecordName = "run.json";
inline constexpr std::string_view kChecksName = "checks.json";

// Where the HP target mass came from; recorded in analysis.json.
enum class HpMassSource { RunRecord, Check };

[[nodiscard]] std::string_view hpMassSourceName(HpMassSource source) noexcept;

struct HpMass {
    double awr = 0.0;             // target mass in neutron masses, as in G4NDL
    double target_mass_MeV = 0.0; // awr * neutron mass
    std::string file;             // G4NDL file the run read
    std::string sha256;
    HpMassSource source = HpMassSource::RunRecord;
    std::string check_id; // empty unless source is Check
};

struct RunInfo {
    std::filesystem::path run_dir;
    std::string run_id;
    int run_record_version = 0;
    std::int32_t schema_version = 0;
    std::uint64_t events = 0;
    std::uint64_t seed = 0;
    std::string git_describe;
    std::vector<std::filesystem::path> event_files; // absolute, in run.json order
    HpMass hp;
};

// Parses an hp_elastic_H1 block. Throws IoError unless awr and target_mass_MeV
// are positive, finite and agree within kHpMassRelTolerance.
[[nodiscard]] HpMass hpMassFromBlock(const nlohmann::json& block, HpMassSource source,
                                     std::string check_id);

// The pure part of readRunInfo: `record` is run.json, `check` the checks.json
// named by `hp_check_id`, if any. With run.json version >= 2 the mass is the
// run's own; a check, when also given, must then report the same block. With
// version 1 the check is required. Throws IoError on any inconsistency.
[[nodiscard]] RunInfo runInfoFromJson(const std::filesystem::path& run_dir,
                                      const nlohmann::json& record,
                                      const std::optional<std::string>& hp_check_id,
                                      const std::optional<nlohmann::json>& check);

// Reads <data_root>/sim/<run_id>/run.json and, with hp_check_id,
// <data_root>/validation/<run_id>/<hp_check_id>/checks.json.
[[nodiscard]] RunInfo readRunInfo(const std::filesystem::path& data_root, std::string_view run_id,
                                  const std::optional<std::string>& hp_check_id);

// Reads a JSON file; throws IoError if it is missing or malformed.
[[nodiscard]] nlohmann::json readJsonFile(const std::filesystem::path& path);

} // namespace riptide::io
