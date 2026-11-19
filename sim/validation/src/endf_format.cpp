#include "riptide/sim/endf_format.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <iterator>
#include <system_error>

namespace riptide::sim::endf {

namespace {

[[noreturn]] void fail(std::string_view what, std::string_view field) {
    throw EndfError(std::string{what} + ": '" + std::string{field} + "'");
}

std::string_view trim(std::string_view text) {
    const auto first = text.find_first_not_of(' ');
    if (first == std::string_view::npos) {
        return {};
    }
    const auto last = text.find_last_not_of(' ');
    return text.substr(first, last - first + 1);
}

bool isExponentLetter(char c) {
    return c == 'e' || c == 'E';
}

} // namespace

double parseEndfFloat(std::string_view field) {
    std::string_view text = trim(field);
    if (text.empty()) {
        return 0.0;
    }
    // std::from_chars accepts a leading '-' but not '+'.
    if (text.front() == '+') {
        text.remove_prefix(1);
    }
    std::string normalized{text};
    if (std::none_of(normalized.begin(), normalized.end(), isExponentLetter)) {
        // A sign after the first character starts the exponent ("1.5-3" is
        // 1.5e-3); the first character may be the mantissa sign.
        const auto sign = normalized.find_first_of("+-", 1);
        if (sign != std::string::npos) {
            normalized.insert(sign, 1, 'e');
        }
    }
    double value = 0.0;
    // from_chars is locale-independent, unlike strtod/stod.
    const char* const begin = normalized.data();
    const char* const end = std::next(begin, static_cast<std::ptrdiff_t>(normalized.size()));
    const auto [ptr, ec] = std::from_chars(begin, end, value);
    if (ec != std::errc{} || ptr != end || !std::isfinite(value)) {
        fail("not an ENDF number", field);
    }
    return value;
}

int parseEndfInt(std::string_view field) {
    const std::string_view text = trim(field);
    if (text.empty()) {
        return 0;
    }
    int value = 0;
    const char* const begin = text.data();
    const char* const end = std::next(begin, static_cast<std::ptrdiff_t>(text.size()));
    const auto [ptr, ec] = std::from_chars(begin, end, value);
    if (ec != std::errc{} || ptr != end) {
        fail("not an ENDF integer", field);
    }
    return value;
}

EndfLine parseEndfLine(std::string_view line) {
    // Tolerate the CR of files written on Windows.
    if (!line.empty() && line.back() == '\r') {
        line.remove_suffix(1);
    }
    if (line.size() < kMinLineWidth) {
        fail("ENDF line shorter than 75 columns", line);
    }
    constexpr std::size_t kMatWidth = 4;
    constexpr std::size_t kMfWidth = 2;
    constexpr std::size_t kMtWidth = 3;
    EndfLine parsed;
    for (std::size_t i = 0; i < kFieldsPerLine; ++i) {
        parsed.fields.at(i) = std::string{line.substr(i * kFieldWidth, kFieldWidth)};
    }
    parsed.mat = parseEndfInt(line.substr(kDataWidth, kMatWidth));
    parsed.mf = parseEndfInt(line.substr(kDataWidth + kMatWidth, kMfWidth));
    parsed.mt = parseEndfInt(line.substr(kDataWidth + kMatWidth + kMfWidth, kMtWidth));
    return parsed;
}

InterpolationLaw interpolationLaw(int code) {
    constexpr int kFirst = static_cast<int>(InterpolationLaw::Histogram);
    constexpr int kLast = static_cast<int>(InterpolationLaw::LogLog);
    if (code < kFirst || code > kLast) {
        throw EndfError("unsupported ENDF interpolation law INT=" + std::to_string(code));
    }
    return static_cast<InterpolationLaw>(code);
}

double interpolate(InterpolationLaw law, double x1, double y1, double x2, double y2, double x) {
    switch (law) {
    case InterpolationLaw::Histogram:
        return y1;
    case InterpolationLaw::LinLin:
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
    case InterpolationLaw::LinLog:
        if (x1 <= 0.0) {
            throw EndfError("lin-log interpolation needs x > 0");
        }
        return y1 + (y2 - y1) * std::log(x / x1) / std::log(x2 / x1);
    case InterpolationLaw::LogLin:
        if (y1 <= 0.0 || y2 <= 0.0) {
            throw EndfError("log-lin interpolation needs y > 0");
        }
        return y1 * std::exp(std::log(y2 / y1) * (x - x1) / (x2 - x1));
    case InterpolationLaw::LogLog:
        if (x1 <= 0.0 || y1 <= 0.0 || y2 <= 0.0) {
            throw EndfError("log-log interpolation needs x > 0 and y > 0");
        }
        return y1 * std::exp(std::log(y2 / y1) * std::log(x / x1) / std::log(x2 / x1));
    }
    throw EndfError("unknown interpolation law");
}

void validateGrid(const std::vector<int>& boundaries, const std::vector<InterpolationLaw>& laws,
                  const std::vector<double>& x, std::string_view where) {
    const auto bad = [where](std::string_view what) {
        throw EndfError(std::string{where} + ": " + std::string{what});
    };
    const std::size_t points = x.size();
    if (points < 2) {
        bad("a table needs at least two points");
    }
    if (boundaries.empty() || boundaries.size() != laws.size()) {
        bad("a table needs one interpolation law per range");
    }
    if (!std::is_sorted(boundaries.begin(), boundaries.end()) ||
        std::adjacent_find(boundaries.begin(), boundaries.end()) != boundaries.end() ||
        boundaries.front() < 2 || static_cast<std::size_t>(boundaries.back()) != points) {
        bad("range boundaries must increase and end at the last point");
    }
    for (std::size_t i = 1; i < points; ++i) {
        if (x[i] < x[i - 1]) {
            bad("x values must not decrease");
        }
        if (i >= 2 && !(x[i] > x[i - 2])) {
            bad("an x value may appear at most twice (a discontinuity)");
        }
    }
}

void validate(const Tab1& table, std::string_view where) {
    validateGrid(table.boundaries, table.laws, table.x, where);
    if (table.y.size() != table.x.size()) {
        throw EndfError(std::string{where} + ": a table needs one y per x");
    }
}

Bracket bracket(const std::vector<int>& boundaries, const std::vector<InterpolationLaw>& laws,
                const std::vector<double>& x, double value) {
    // Written as a negated range test so that NaN is refused too.
    // NOLINTNEXTLINE(readability-simplify-boolean-expr)
    if (!(value >= x.front() && value <= x.back())) {
        throw EndfError("x = " + std::to_string(value) + " outside the tabulated range [" +
                        std::to_string(x.front()) + ", " + std::to_string(x.back()) + "]");
    }
    // First point strictly above the value: the interval below it starts at
    // the last point <= value, which skips the zero-width interval of a
    // discontinuity.
    const auto above = std::upper_bound(x.begin(), x.end(), value);
    if (above == x.end()) {
        const std::size_t last = x.size() - 1;
        return {last, last, laws.back()};
    }
    const auto upper = static_cast<std::size_t>(std::distance(x.begin(), above));
    // The interval (upper - 1, upper) belongs to the first range whose 1-based
    // boundary reaches the upper point.
    const auto range =
        std::lower_bound(boundaries.begin(), boundaries.end(), static_cast<int>(upper + 1));
    const auto law = laws.at(static_cast<std::size_t>(std::distance(boundaries.begin(), range)));
    return {upper - 1, upper, law};
}

double evaluate(const Tab1& table, double x) {
    const Bracket where = bracket(table.boundaries, table.laws, table.x, x);
    if (where.lower == where.upper) {
        return table.y[where.lower];
    }
    return interpolate(where.law, table.x[where.lower], table.y[where.lower], table.x[where.upper],
                       table.y[where.upper], x);
}

} // namespace riptide::sim::endf
