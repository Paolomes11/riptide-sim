#pragma once

// Where simulation and analysis data live. Data never go inside the repository
// nor, by default, on the system disk: the root is an external SSD, chosen with
// precedence command line > RIPTIDE_DATA_DIR > default. Nothing here falls back
// to another location; a missing or unsuitable root is an error the user must
// resolve (mount the disk, or pass an explicit directory).
//
// Layout below the root:
//   sim/<run_id>/                     events_t*.root, run.json
//   analysis/<run_id>/<analysis_id>/  results, analysis.json
//   validation/<run_id>/<check_id>/   physics checks of a run, checks.json

#include <sys/stat.h>
#include <unistd.h>

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace riptide::data {

inline constexpr std::string_view kDataDirEnvVar = "RIPTIDE_DATA_DIR";
inline constexpr std::string_view kDefaultDataRoot = "/mnt/external_ssd/riptide-sim";
// Headroom left free on the target disk so a run never fills it completely.
inline constexpr std::uintmax_t kFreeSpaceMargin = std::uintmax_t{2} << 30U; // 2 GiB
inline constexpr std::size_t kMaxIdLength = 128;

class DataPathError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Resolves the data root. `env_value` is the content of RIPTIDE_DATA_DIR (or
// nullptr when unset), passed in so the policy is testable without touching the
// process environment. Empty strings count as unset. Returns an absolute path.
inline std::filesystem::path resolveDataRoot(const std::optional<std::string>& cli_value,
                                             const char* env_value) {
    std::filesystem::path root;
    if (cli_value && !cli_value->empty()) {
        root = *cli_value;
    } else if (env_value != nullptr && *env_value != '\0') {
        root = env_value;
    } else {
        root = std::string{kDefaultDataRoot};
    }
    // A relative root would silently depend on the working directory.
    if (!root.is_absolute()) {
        throw DataPathError("data root must be an absolute path: " + root.string());
    }
    return root.lexically_normal();
}

// Same as above, reading RIPTIDE_DATA_DIR from the environment.
inline std::filesystem::path resolveDataRoot(const std::optional<std::string>& cli_value) {
    // NOLINTNEXTLINE(concurrency-mt-unsafe): called once at startup, before any thread
    return resolveDataRoot(cli_value, std::getenv(std::string{kDataDirEnvVar}.c_str()));
}

// True when both existing paths are on the same filesystem (same st_dev).
inline bool sameDevice(const std::filesystem::path& a, const std::filesystem::path& b) {
    struct stat sa {};
    struct stat sb {};
    if (::stat(a.c_str(), &sa) != 0) {
        throw DataPathError("cannot stat " + a.string());
    }
    if (::stat(b.c_str(), &sb) != 0) {
        throw DataPathError("cannot stat " + b.string());
    }
    return sa.st_dev == sb.st_dev;
}

// Bytes needed by a run, with an explicit error instead of a silent overflow.
inline std::uintmax_t estimateRunBytes(std::uint64_t n_events, std::uint64_t bytes_per_event) {
    if (bytes_per_event != 0 &&
        n_events > std::numeric_limits<std::uintmax_t>::max() / bytes_per_event) {
        throw DataPathError("run size estimate overflows");
    }
    return std::uintmax_t{n_events} * bytes_per_event;
}

// Checks that `root` can receive `required_bytes` more data: it must exist, be a
// writable directory, not be on the system disk (unless explicitly allowed) and
// have at least required_bytes + margin_bytes available. Throws DataPathError
// with the reason otherwise. The root is never created here: if the SSD is not
// mounted its directory is missing, and creating it would write to the system
// disk under the mount point.
inline void requireWritableRoot(const std::filesystem::path& root, std::uintmax_t required_bytes,
                                std::uintmax_t margin_bytes = kFreeSpaceMargin,
                                bool allow_system_disk = false) {
    std::error_code ec;
    if (!std::filesystem::is_directory(root, ec)) {
        throw DataPathError(
            "data root does not exist or is not a directory (is the SSD mounted?): " +
            root.string());
    }
    if (::access(root.c_str(), W_OK) != 0) {
        throw DataPathError("data root is not writable: " + root.string());
    }
    if (!allow_system_disk && sameDevice(root, "/")) {
        throw DataPathError("data root is on the system disk: " + root.string() +
                            " (mount the external SSD or allow it explicitly)");
    }
    const std::filesystem::space_info info = std::filesystem::space(root, ec);
    if (ec) {
        throw DataPathError("cannot query free space of " + root.string() + ": " + ec.message());
    }
    if (required_bytes > std::numeric_limits<std::uintmax_t>::max() - margin_bytes ||
        info.available < required_bytes + margin_bytes) {
        throw DataPathError("not enough free space in " + root.string() + ": available " +
                            std::to_string(info.available) + " B, required " +
                            std::to_string(required_bytes) + " B + margin " +
                            std::to_string(margin_bytes) + " B");
    }
}

// Run and analysis identifiers become directory names: only [A-Za-z0-9._-],
// not starting with '.' or '-', so they cannot escape the root or look like options.
inline bool isValidId(std::string_view id) {
    if (id.empty() || id.size() > kMaxIdLength || id.front() == '.' || id.front() == '-') {
        return false;
    }
    for (const char c : id) {
        const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                        (c >= '0' && c <= '9') || c == '.' || c == '_' || c == '-';
        if (!ok) {
            return false;
        }
    }
    return true;
}

inline void requireValidId(std::string_view id, std::string_view what) {
    if (!isValidId(id)) {
        throw DataPathError("invalid " + std::string{what} + ": '" + std::string{id} + "'");
    }
}

// <root>/sim/<run_id>
inline std::filesystem::path simRunDir(const std::filesystem::path& root, std::string_view run_id) {
    requireValidId(run_id, "run id");
    return root / "sim" / std::string{run_id};
}

// <root>/analysis/<run_id>/<analysis_id>
inline std::filesystem::path analysisDir(const std::filesystem::path& root, std::string_view run_id,
                                         std::string_view analysis_id) {
    requireValidId(run_id, "run id");
    requireValidId(analysis_id, "analysis id");
    return root / "analysis" / std::string{run_id} / std::string{analysis_id};
}

// <root>/validation/<run_id>/<check_id>
inline std::filesystem::path validationDir(const std::filesystem::path& root,
                                           std::string_view run_id, std::string_view check_id) {
    requireValidId(run_id, "run id");
    requireValidId(check_id, "check id");
    return root / "validation" / std::string{run_id} / std::string{check_id};
}

// <root>/stats/<calibration_id>: outputs of the statistics stage that do not
// belong to one simulation run (pseudo-experiment calibrations).
inline std::filesystem::path statsDir(const std::filesystem::path& root,
                                      std::string_view calibration_id) {
    requireValidId(calibration_id, "calibration id");
    return root / "stats" / std::string{calibration_id};
}

} // namespace riptide::data
