#pragma once

// Pieces shared by the command lines of the riptide_sim applications.

#include <charconv>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace riptide::sim {

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
        throw CliError(std::string{option} + ": expected a non-negative integer, got '" + text +
                       "'");
    }
    return value;
}

} // namespace riptide::sim
