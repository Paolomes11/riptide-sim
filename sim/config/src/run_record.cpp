#include "riptide/sim/run_record.hpp"

#include <nlohmann/json.hpp>

#include <array>
#include <fstream>
#include <system_error>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/version.hpp"

namespace riptide::sim {

namespace {

constexpr int kRunRecordIndent = 2;

// "1970-01-01T00:00:00Z" plus the terminating null.
constexpr std::size_t kUtcBufferSize = 21;

// Suffix of the temporary file renamed onto run.json.
constexpr const char* kTemporarySuffix = ".tmp";

} // namespace

std::string formatUtc(std::time_t time) {
    std::tm utc{};
    // gmtime_r, not gmtime: the latter returns a pointer to shared static storage.
    if (gmtime_r(&time, &utc) == nullptr) {
        throw RunRecordError("cannot convert time to UTC");
    }
    std::array<char, kUtcBufferSize> buffer{};
    if (std::strftime(buffer.data(), buffer.size(), "%Y-%m-%dT%H:%M:%SZ", &utc) == 0) {
        throw RunRecordError("cannot format UTC time");
    }
    return std::string{buffer.data()};
}

std::vector<std::string> eventFileNames(int threads) {
    // Mirrors the Geant4 rule for unmerged ntuples: "<base>_t<thread id>.<ext>".
    const std::filesystem::path base{kEventFileName};
    const auto stem = base.stem().string();
    const auto extension = base.extension().string();
    std::vector<std::string> names;
    for (int k = 0; k < threads; ++k) {
        names.push_back(stem + "_t" + std::to_string(k) + extension);
    }
    return names;
}

nlohmann::json runRecordJson(const RunConfig& config, const RunProvenance& run) {
    nlohmann::json doc;
    doc["run_record_version"] = kRunRecordVersion;
    doc["run_id"] = run.run_id;
    doc["project_version"] = std::string{kProjectVersion};
    doc["git_describe"] = std::string{kGitDescribe};
    doc["geant4_version"] = run.geant4_version;
    doc["schema_version"] = schema::kSchemaVersion;
    doc["config_version"] = config.config_version;
    doc["seed"] = run.seed;
    doc["events"] = run.events;
    doc["threads"] = run.threads;
    doc["started_utc"] = run.started_utc;
    doc["elapsed_s"] = run.elapsed_s;
    doc["event_files"] = eventFileNames(run.threads);
    doc["config"] = toJson(config);
    return doc;
}

void writeRunRecord(const std::filesystem::path& run_dir, const nlohmann::json& doc) {
    const auto target = run_dir / kRunRecordFileName;
    if (std::filesystem::exists(target)) {
        throw RunRecordError("run record already exists: " + target.string());
    }
    auto temporary = target;
    temporary += kTemporarySuffix;
    {
        std::ofstream out(temporary, std::ios::out | std::ios::trunc);
        out << doc.dump(kRunRecordIndent) << '\n';
        out.close();
        if (!out) {
            throw RunRecordError("cannot write " + temporary.string());
        }
    }
    // rename is atomic within one filesystem: readers see either no run.json or
    // the complete file, never a partial one.
    std::error_code error;
    std::filesystem::rename(temporary, target, error);
    if (error) {
        throw RunRecordError("cannot rename " + temporary.string() + " to " + target.string() +
                             ": " + error.message());
    }
}

} // namespace riptide::sim
