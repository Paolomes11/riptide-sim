#pragma once

// analysis.json and cutflow.csv: what an analysis read, how it selected and
// smeared the events, and how many survived each cut. The record has no
// timestamps, so the same run, settings and code give a byte-identical file.

#include <nlohmann/json.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

#include "riptide/io/response.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"

namespace riptide::io {

// Layout of analysis.json, incremented when its keys change.
inline constexpr int kAnalysisRecordVersion = 1;
inline constexpr std::string_view kAnalysisRecordName = "analysis.json";
inline constexpr std::string_view kCutFlowName = "cutflow.csv";

struct AnalysisSettings {
    std::string analysis_id;
    std::uint64_t seed = 0;
    Resolution resolution;
    Selection selection;
    std::uint64_t max_events = 0; // 0 = all events of the run
};

// `summaries` holds the residual summaries computed by the application; it is
// stored as given under "summaries".
[[nodiscard]] nlohmann::json analysisRecordJson(const RunInfo& run,
                                                const AnalysisSettings& settings,
                                                const CutFlow& cut_flow,
                                                const nlohmann::json& summaries);

// One line per cut: index, name, threshold, events passing up to that cut,
// fraction of the total. The first line counts all events read.
[[nodiscard]] std::string cutFlowCsv(const Selection& selection, const CutFlow& cut_flow);

// Writes `content` to `target` through a temporary file and a rename, so that
// a reader sees either no file or the complete one. Refuses to overwrite.
void writeFileAtomically(const std::filesystem::path& target, const std::string& content);

} // namespace riptide::io
