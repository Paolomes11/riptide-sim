#pragma once

// run.json: the provenance record of one simulation run. The Geant4 ntuple has
// no attributes, so everything a reader needs besides the rows lives here: the
// effective configuration, the seed, the software versions and the schema
// version of the event files.
//
// run.json is also the completion marker of a run. It is written atomically
// (temporary file, then rename) only after BeamOn has returned, so a run
// directory without it is an interrupted run and readers refuse it.

#include <nlohmann/json_fwd.hpp>

#include <cstdint>
#include <ctime>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "riptide/sim/run_config.hpp"

namespace riptide::sim {

// Version of the run.json layout; bumped on any change of its keys or meaning.
inline constexpr int kRunRecordVersion = 1;

inline constexpr const char* kRunRecordFileName = "run.json";

// Base name handed to the Geant4 analysis manager. With ntuple merging off,
// each worker writes "events_t<thread id>.root" and the master writes nothing.
inline constexpr const char* kEventFileName = "events.root";

class RunRecordError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct RunProvenance {
    std::string run_id;
    std::uint64_t seed = 0;
    std::uint64_t events = 0;
    int threads = 0;
    std::string geant4_version; // "major.minor.patch"
    std::string started_utc;    // ISO 8601, see formatUtc
    double elapsed_s = 0.0;     // wall-clock time from Initialize to the end of BeamOn [s]
};

// "YYYY-MM-DDTHH:MM:SSZ"; throws RunRecordError if the time cannot be converted.
[[nodiscard]] std::string formatUtc(std::time_t time);

// Per-thread event files of a run with `threads` workers, in thread order.
[[nodiscard]] std::vector<std::string> eventFileNames(int threads);

// The full run.json document; the configuration is stored with toJson, so
// parseRunConfig(doc["config"]) gives back `config`.
[[nodiscard]] nlohmann::json runRecordJson(const RunConfig& config, const RunProvenance& run);

// Writes run_dir/run.json atomically. Throws RunRecordError if run.json already
// exists (a completed run is never overwritten) or if writing fails.
void writeRunRecord(const std::filesystem::path& run_dir, const nlohmann::json& doc);

} // namespace riptide::sim
