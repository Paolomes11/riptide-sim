#include "riptide/io/selection.hpp"

#include <charconv>
#include <cmath>
#include <set>
#include <sstream>
#include <system_error>

#include "riptide/common/event_schema.hpp"
#include "riptide/io/run_info.hpp"

namespace riptide::io {

namespace {

constexpr std::string_view kSingleNpName = "single_np";
constexpr std::string_view kContainedName = "contained";
constexpr std::string_view kEvisMinName = "evis_min";
constexpr char kCutSeparator = ',';
constexpr char kValueSeparator = '=';
// Digits enough to print any double so that it parses back to the same value.
constexpr int kRoundTripDigits = 17;

double parseThreshold(std::string_view text) {
    double value = 0.0;
    const auto* const end = text.data() + text.size();
    const auto [last, error] = std::from_chars(text.data(), end, value);
    if (text.empty() || error != std::errc{} || last != end) {
        throw IoError("evis_min: malformed threshold '" + std::string{text} + "'");
    }
    if (!std::isfinite(value) || value < 0.0) {
        throw IoError("evis_min: threshold must be finite and >= 0 MeV");
    }
    return value;
}

CutSpec parseCut(std::string_view item) {
    const auto equals = item.find(kValueSeparator);
    const auto name = item.substr(0, equals);
    const bool has_value = equals != std::string_view::npos;
    if (name == kEvisMinName) {
        if (!has_value) {
            throw IoError("evis_min needs an explicit threshold: evis_min=<MeV>");
        }
        return {Cut::EvisMin, parseThreshold(item.substr(equals + 1))};
    }
    if (has_value) {
        throw IoError("cut '" + std::string{name} + "' takes no value");
    }
    if (name == kSingleNpName) {
        return {Cut::SingleNp, 0.0};
    }
    if (name == kContainedName) {
        return {Cut::Contained, 0.0};
    }
    throw IoError("unknown cut '" + std::string{item} + "'");
}

bool passes(const CutSpec& spec, const CutInputs& event) {
    switch (spec.cut) {
    case Cut::SingleNp:
        return event.truth_event_class == static_cast<std::int32_t>(schema::EventClass::SingleNp);
    case Cut::Contained:
        return event.truth_p_exits == 0;
    case Cut::EvisMin:
        return event.det_Evis > spec.threshold;
    }
    return false;
}

} // namespace

std::string_view cutName(Cut cut) noexcept {
    switch (cut) {
    case Cut::SingleNp:
        return kSingleNpName;
    case Cut::Contained:
        return kContainedName;
    case Cut::EvisMin:
        return kEvisMinName;
    }
    return "unknown";
}

Selection Selection::parse(std::string_view text) {
    Selection selection;
    std::set<Cut> seen;
    std::size_t start = 0;
    while (start <= text.size()) {
        const auto comma = text.find(kCutSeparator, start);
        const auto stop = comma == std::string_view::npos ? text.size() : comma;
        const auto item = text.substr(start, stop - start);
        if (item.empty()) {
            throw IoError("empty cut in selection '" + std::string{text} + "'");
        }
        const auto spec = parseCut(item);
        if (!seen.insert(spec.cut).second) {
            throw IoError("cut '" + std::string{cutName(spec.cut)} + "' appears twice");
        }
        selection.cuts_.push_back(spec);
        start = stop + 1;
    }
    return selection;
}

std::optional<std::size_t> Selection::firstFailingCut(const CutInputs& event) const {
    for (std::size_t k = 0; k < cuts_.size(); ++k) {
        if (!passes(cuts_[k], event)) {
            return k;
        }
    }
    return std::nullopt;
}

std::string Selection::toString() const {
    std::ostringstream out;
    out.precision(kRoundTripDigits);
    for (std::size_t k = 0; k < cuts_.size(); ++k) {
        if (k > 0) {
            out << kCutSeparator;
        }
        out << cutName(cuts_[k].cut);
        if (cuts_[k].cut == Cut::EvisMin) {
            out << kValueSeparator << cuts_[k].threshold;
        }
    }
    return out.str();
}

void CutFlow::record(std::optional<std::size_t> first_failing) {
    ++total_;
    const auto reached = first_failing.value_or(passed_.size());
    if (reached > passed_.size()) {
        throw IoError("cut index out of range in the cut flow");
    }
    for (std::size_t k = 0; k < reached; ++k) {
        ++passed_[k];
    }
}

} // namespace riptide::io
