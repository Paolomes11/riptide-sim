#include "riptide/io/acceptance_record.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <string>
#include <utility>

#include "riptide/io/analysis_record.hpp"

namespace riptide::io {

namespace {

// A last bin shorter than this fraction of the step is rounding, not a bin.
constexpr double kSliverFraction = 1e-9;

bool isFinite(double value) {
    return std::isfinite(value);
}

void requireRange(double lo, double hi) {
    if (!isFinite(lo) || !isFinite(hi) || !(lo < hi)) {
        throw IoError("edges need finite lo < hi");
    }
}

void requireEdges(const std::vector<double>& edges, const char* what) {
    if (edges.size() < 2) {
        throw IoError(std::string{what} + " edges need at least one bin");
    }
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (!isFinite(edges[i]) || (i > 0 && !(edges[i - 1] < edges[i]))) {
            throw IoError(std::string{what} + " edges must be finite and strictly increasing");
        }
    }
}

// Bin of `value`: [e_i, e_{i+1}) except the last, closed at the top.
std::size_t binOf(const std::vector<double>& edges, double value, const char* what) {
    if (!(value >= edges.front() && value <= edges.back())) {
        throw IoError(std::string{"an event has "} + what + " outside the acceptance edges");
    }
    const auto upper = std::upper_bound(edges.begin(), edges.end(), value);
    const auto bin = static_cast<std::size_t>(std::distance(edges.begin(), upper)) - 1;
    return std::min(bin, edges.size() - 2);
}

bool inHalf(std::int32_t event_id, EventHalf half) {
    switch (half) {
    case EventHalf::Even:
        return event_id % 2 == 0;
    case EventHalf::Odd:
        return event_id % 2 != 0;
    case EventHalf::All:
        break;
    }
    return true;
}

} // namespace

std::string_view eventHalfName(EventHalf half) noexcept {
    switch (half) {
    case EventHalf::Even:
        return "even";
    case EventHalf::Odd:
        return "odd";
    case EventHalf::All:
        break;
    }
    return "all";
}

EventHalf parseEventHalf(std::string_view text) {
    for (const EventHalf half : {EventHalf::All, EventHalf::Even, EventHalf::Odd}) {
        if (text == eventHalfName(half)) {
            return half;
        }
    }
    throw IoError("event half must be all, even or odd, not '" + std::string{text} + "'");
}

std::vector<double> steppedEdges(double lo, double hi, double step) {
    requireRange(lo, hi);
    if (!isFinite(step) || !(step > 0.0)) {
        throw IoError("the edge step must be positive and finite");
    }
    std::vector<double> edges{lo};
    // Edges from the integer multiple, not by repeated addition, so rounding
    // does not accumulate.
    for (std::size_t i = 1;; ++i) {
        const double edge = lo + (static_cast<double>(i) * step);
        if (edge >= hi - (kSliverFraction * step)) {
            break;
        }
        edges.push_back(edge);
    }
    edges.push_back(hi);
    return edges;
}

std::vector<double> uniformEdges(double lo, double hi, std::size_t bins) {
    requireRange(lo, hi);
    if (bins == 0) {
        throw IoError("uniform edges need at least one bin");
    }
    std::vector<double> edges(bins + 1);
    const auto n = static_cast<double>(bins);
    for (std::size_t i = 0; i < bins; ++i) {
        edges[i] = lo + ((hi - lo) * static_cast<double>(i) / n);
    }
    edges[bins] = hi;
    return edges;
}

AcceptanceTable fillAcceptanceTable(const EventSample& sample, std::vector<double> En_edges,
                                    std::vector<double> theta_edges, EventHalf half) {
    requireEdges(En_edges, "E_n");
    requireEdges(theta_edges, "theta_p");
    const std::size_t n = sample.size();
    if (sample.truth_En.size() != n || sample.truth_theta_p.size() != n) {
        throw IoError("truth_En and truth_theta_p need one entry per event");
    }
    AcceptanceTable table;
    const std::size_t n_theta = theta_edges.size() - 1;
    table.counts.assign((En_edges.size() - 1) * n_theta, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (!inHalf(sample.event_id[i], half)) {
            continue;
        }
        const std::size_t k = binOf(En_edges, sample.truth_En[i], "truth_En");
        const std::size_t t = binOf(theta_edges, sample.truth_theta_p[i], "truth_theta_p");
        table.counts[(k * n_theta) + t] += 1.0;
    }
    table.En_edges = std::move(En_edges);
    table.theta_edges = std::move(theta_edges);
    return table;
}

nlohmann::json acceptanceRecordJson(std::string_view stats_id, const RunInfo& run,
                                    const Selection& selection, const CutFlow& cut_flow,
                                    EventHalf half, const AcceptanceTable& table) {
    nlohmann::json doc;
    doc["acceptance_version"] = kAcceptanceRecordVersion;
    doc["stats_id"] = std::string{stats_id};
    doc["source"] = {{"run_id", run.run_id},
                     {"run_seed", run.seed},
                     {"run_record_version", run.run_record_version},
                     {"schema_version", run.schema_version},
                     {"run_git_describe", run.git_describe}};
    doc["cuts"] = selection.toString();
    doc["cut_flow"] = cutsJson(selection, cut_flow);
    doc["events_read"] = cut_flow.total();
    doc["events_selected"] = cut_flow.selected();
    doc["half"] = std::string{eventHalfName(half)};
    double entries = 0.0;
    for (const double c : table.counts) {
        entries += c;
    }
    doc["entries"] = entries;
    doc["En_edges_MeV"] = table.En_edges;
    doc["theta_edges_rad"] = table.theta_edges;
    doc["counts"] = table.counts;
    return doc;
}

AcceptanceTable acceptanceTableFromRecord(const nlohmann::json& record) {
    try {
        if (record.at("acceptance_version").get<int>() != kAcceptanceRecordVersion) {
            throw IoError("acceptance.json has version " + record.at("acceptance_version").dump() +
                          ", expected " + std::to_string(kAcceptanceRecordVersion));
        }
        AcceptanceTable table;
        table.En_edges = record.at("En_edges_MeV").get<std::vector<double>>();
        table.theta_edges = record.at("theta_edges_rad").get<std::vector<double>>();
        table.counts = record.at("counts").get<std::vector<double>>();
        requireEdges(table.En_edges, "E_n");
        requireEdges(table.theta_edges, "theta_p");
        if (table.counts.size() != (table.En_edges.size() - 1) * (table.theta_edges.size() - 1)) {
            throw IoError("acceptance.json counts do not match its edges");
        }
        return table;
    } catch (const nlohmann::json::exception& e) {
        throw IoError(std::string{"acceptance.json: "} + e.what());
    }
}

} // namespace riptide::io
