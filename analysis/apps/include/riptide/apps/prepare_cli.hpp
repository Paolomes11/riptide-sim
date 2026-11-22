#pragma once

// Command line of riptide_analysis_prepare. Parsing is a pure function of the
// argument list, unit-tested in the CI build; the application then resolves
// the data root and checks the directories.
//
// Every analysis setting is required: the seed, the resolutions and the cuts
// define the analysis, and a silent default would make two analyses look
// alike in their records while differing in what they did.

#include <charconv>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "riptide/io/response.hpp"
#include "riptide/io/selection.hpp"

namespace riptide::apps {

class CliError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Numbers are taken as strings and converted here: stream extraction into an
// unsigned type accepts "-1" and wraps it to a huge value, and std::from_chars
// rejects signs, trailing characters and overflow. Throws CliError naming
// `option` on any malformed value.
template <typename T>
[[nodiscard]] T parseCliNumber(const std::string& text, std::string_view option) {
    T value{};
    const char* first = text.data();
    const char* last = first + text.size();
    const auto [ptr, ec] = std::from_chars(first, last, value);
    if (text.empty() || ec != std::errc{} || ptr != last) {
        throw CliError(std::string{option} + ": malformed number '" + text + "'");
    }
    return value;
}

struct PrepareOptions {
    bool show_help = false;
    bool show_version = false;
    std::optional<std::string> data_dir; // unset: RIPTIDE_DATA_DIR, then the default root
    std::string run_id;
    std::string analysis_id;
    std::uint64_t seed = 0;
    io::Resolution resolution;
    io::Selection selection;
    std::optional<std::string> hp_check_id;
    std::uint64_t max_events = 0; // 0 = all events of the run
};

// `args` excludes the program name. Throws CliError on unknown options,
// malformed numbers, negative or non-finite resolutions, invalid ids or cuts
// and missing required options. With --help or --version the other options
// are not required.
[[nodiscard]] PrepareOptions parsePrepareCommandLine(const std::vector<std::string>& args);

// Usage text printed by --help.
[[nodiscard]] std::string prepareHelpText();

} // namespace riptide::apps
