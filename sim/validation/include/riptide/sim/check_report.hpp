#pragma once

// Output of riptide_sim_check: the results as JSON (the summary that goes into
// checks.json) and the per-bin tables as CSV. Pure formatting, no file system
// access except the atomic write, so the layout is unit-tested in every preset.
// Observed and expected values are always written side by side and unmerged;
// the merged bins appear only in the JSON, next to the chi^2 they produced.

#include <nlohmann/json.hpp>

#include <cstdint>
#include <filesystem>
#include <ostream>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "riptide/sim/check_accumulator.hpp"

namespace riptide::sim::checks {

inline constexpr const char* kChecksFileName = "checks.json";
inline constexpr const char* kSpectrumUFileName = "spectrum_u.csv";
inline constexpr const char* kDepthZFileName = "depth_z.csv";
inline constexpr const char* kFractionsFileName = "fractions.csv";
inline constexpr const char* kLegendreFileName = "legendre.csv";
inline constexpr const char* kResidualsFileName = "residuals_Ep.csv";

// Version of the checks.json layout; bumped on any change of its keys.
inline constexpr int kChecksVersion = 1;

class ReportError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Name of a FirstInteraction code in the reports, e.g. "np_elastic".
// Throws std::out_of_range for a code outside the schema.
[[nodiscard]] std::string_view categoryName(std::int32_t code);

// The "results" object of checks.json: event counts, the five checks per
// library and the consistency counters, the must-be-zero ones apart from the
// expected-rare ones. A skipped check carries its reason in "skipped" and no
// chi^2; a compared one has "skipped": null.
[[nodiscard]] nlohmann::json resultsJson(const CheckResults& results);

// low,high,observed,expected_<library>... one row per original bin.
void writeHistogramCsv(std::ostream& out, const HistogramResult& histogram);
// code,category,observed,expected_<library>...
void writeFractionsCsv(std::ostream& out, const FractionResult& fractions);
// l,value,error,expected_<library>,pull_<library>...; only the header when skipped.
void writeLegendreCsv(std::ostream& out, const LegendreResult& legendre);
// mu,Ep_residual_MeV, one row per n-p event in the order the events were read.
// Throws std::invalid_argument if the two columns differ in length.
void writeResidualsCsv(std::ostream& out, const std::vector<double>& mu,
                       const std::vector<double>& Ep_residuals);

// Writes `doc` to `path` through a temporary file and a rename, so a reader
// sees either no file or the complete one: the presence of checks.json marks a
// finished check. Throws ReportError if `path` exists or writing fails.
void writeJsonAtomically(const std::filesystem::path& path, const nlohmann::json& doc);

} // namespace riptide::sim::checks
