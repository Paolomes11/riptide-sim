#pragma once

// acceptance.json (D40): the (E_n, theta_p) histogram of the selected events
// of one flat Geant4 run, from which the Case A model takes its theta_p prior
// per E_n row. The record keeps the counts, not the ratio: the ratio is a
// choice of the stats library, the counts are what the run measured.

#include <nlohmann/json.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "riptide/common/event_sample.hpp"
#include "riptide/common/units.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"

namespace riptide::io {

// Layout of acceptance.json, incremented when its keys change.
inline constexpr int kAcceptanceRecordVersion = 1;
inline constexpr std::string_view kAcceptanceRecordName = "acceptance.json";

// Binning of D40: 1 MeV in E_n (the last bin may be shorter) and 5 degrees in
// theta_p on [0, pi/2].
inline constexpr double kAcceptanceEnStep = 1.0 * units::MeV;
inline constexpr std::size_t kAcceptanceThetaBins = 18;

// Which events of the run enter: all, or half of them by event_id parity, to
// measure how much the table moves with its own statistics.
enum class EventHalf { All, Even, Odd };

[[nodiscard]] std::string_view eventHalfName(EventHalf half) noexcept;

// "all", "even" or "odd"; throws IoError otherwise.
[[nodiscard]] EventHalf parseEventHalf(std::string_view text);

struct AcceptanceTable {
    std::vector<double> En_edges;    // MeV
    std::vector<double> theta_edges; // rad
    std::vector<double> counts;      // [k * n_theta + t]
};

// Edges from lo to hi in steps of `step`, the last bin ending at hi. A last
// bin shorter than a part in 1e9 of the step is merged into the previous one,
// so rounding never leaves a sliver. Throws IoError unless lo < hi are finite
// and step is positive and finite.
[[nodiscard]] std::vector<double> steppedEdges(double lo, double hi, double step);

// `bins` equal bins on [lo, hi], the last edge exactly hi. Throws IoError
// unless lo < hi are finite and bins > 0.
[[nodiscard]] std::vector<double> uniformEdges(double lo, double hi, std::size_t bins);

// Histogram of truth_En and truth_theta_p of the events of `sample` in `half`.
// Bins are [lo, hi) except the last, which is closed so the upper edge of the
// source belongs to it. Throws IoError if an event falls outside the edges,
// if the edges are not strictly increasing with at least one bin, or if the
// columns are not filled.
[[nodiscard]] AcceptanceTable fillAcceptanceTable(const EventSample& sample,
                                                  std::vector<double> En_edges,
                                                  std::vector<double> theta_edges, EventHalf half);

// The record: version, ids, the run it comes from, cuts and cut-flow, the half,
// edges and counts.
[[nodiscard]] nlohmann::json acceptanceRecordJson(std::string_view stats_id, const RunInfo& run,
                                                  const Selection& selection,
                                                  const CutFlow& cut_flow, EventHalf half,
                                                  const AcceptanceTable& table);

// Reads the table of a record written by acceptanceRecordJson. Throws IoError
// on another version, a missing key, or counts that do not match the edges.
[[nodiscard]] AcceptanceTable acceptanceTableFromRecord(const nlohmann::json& record);

} // namespace riptide::io
