#pragma once

// Reading back a riptide_sim run: run.json (provenance and completion marker)
// and the per-thread event files. Rows can be streamed one at a time, so large
// runs never sit in memory; loadRun keeps a whole run only for small ones.
// The implementation uses the Geant4 ROOT reader, but this header does not
// expose Geant4: the simulation side keeps no ROOT dependency.

#include <nlohmann/json.hpp>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <stdexcept>
#include <string_view>

#include "riptide/sim/event_row.hpp"

namespace riptide::sim {

class RunReadError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Reads run_dir/run.json and checks that a reader can trust it: the file
// exists (without it the run was interrupted), run_record_version is the one
// this build writes, schema_version is one it can read (kOldestReadable...
// up to kSchemaVersion), and event_files is a list of names. Throws
// RunReadError otherwise.
[[nodiscard]] nlohmann::json readRunRecord(const std::filesystem::path& run_dir);

// Calls `visit` on every row of every event file listed in `record`, in file
// order and, within a file, in write order. Every column of the record's
// schema version is bound by name, so a file missing one throws RunReadError;
// columns added by later versions are 0 in the rows of an older file. The files are closed on
// return and on exception. Returns the number of rows visited.
std::size_t forEachEventRow(const std::filesystem::path& run_dir, const nlohmann::json& record,
                            const std::function<void(const EventRow&)>& visit);

// A whole run in memory, keyed by event_id. A row costs a few hundred bytes
// here, so this is meant for short runs (tests), not production runs.
struct RunData {
    nlohmann::json record;
    std::map<std::int32_t, EventRow> rows;
    std::size_t row_count = 0; // rows read, duplicated event ids included
};

[[nodiscard]] RunData loadRun(const std::filesystem::path& run_dir);

// Typed access to a column by schema name. An unknown name throws
// std::out_of_range, a column of the other type std::bad_variant_access.
[[nodiscard]] std::int32_t intAt(const EventRow& row, std::string_view name);
[[nodiscard]] double doubleAt(const EventRow& row, std::string_view name);

} // namespace riptide::sim
