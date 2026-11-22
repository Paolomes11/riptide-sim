#include "riptide/io/event_columns.hpp"

#include <algorithm>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

#include "riptide/io/run_info.hpp"

namespace riptide::io {

namespace {

template <typename T>
std::vector<T> permuted(const std::vector<T>& column, const std::vector<std::size_t>& order) {
    std::vector<T> result;
    result.reserve(order.size());
    for (const auto index : order) {
        result.push_back(column[index]);
    }
    return result;
}

} // namespace

std::uint64_t expectedEventCount(std::uint64_t events_in_run, std::uint64_t max_events) noexcept {
    return max_events == 0 ? events_in_run : std::min(events_in_run, max_events);
}

void sortByEventId(EventColumns& columns, std::uint64_t expected) {
    const auto n = columns.size();
    const auto check_length = [n](std::string_view name, const auto& column) {
        if (column.size() != n) {
            throw IoError("column " + std::string{name} + " has " + std::to_string(column.size()) +
                          " entries, event_id has " + std::to_string(n));
        }
    };
    columns.forEachInt32Column(check_length);
    columns.forEachDoubleColumn(check_length);
    if (n != expected) {
        throw IoError("read " + std::to_string(n) + " events, expected " +
                      std::to_string(expected));
    }

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    // Ids are checked to be unique below, so a plain sort is deterministic.
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return columns.event_id[a] < columns.event_id[b];
    });
    for (std::size_t i = 0; i < n; ++i) {
        const auto id = columns.event_id[order[i]];
        if (id < 0 || static_cast<std::size_t>(id) != i) {
            throw IoError("event ids are not 0.." + std::to_string(n) + "-1 each once: position " +
                          std::to_string(i) + " holds id " + std::to_string(id));
        }
    }
    columns.forEachInt32Column([&order](std::string_view, std::vector<std::int32_t>& column) {
        column = permuted(column, order);
    });
    columns.forEachDoubleColumn([&order](std::string_view, std::vector<double>& column) {
        column = permuted(column, order);
    });
}

EventSample selectEvents(const EventColumns& columns, const Selection& selection,
                         CutFlow& cut_flow) {
    if (cut_flow.cuts() != selection.cuts().size()) {
        throw IoError("the cut flow does not belong to the selection");
    }
    EventSample sample;
    for (std::size_t i = 0; i < columns.size(); ++i) {
        const auto failing =
            selection.firstFailingCut({columns.truth_event_class[i], columns.truth_p_exits[i],
                                       columns.det_Evis[i], columns.truth_valid_np[i]});
        cut_flow.record(failing);
        if (failing.has_value()) {
            continue;
        }
        sample.event_id.push_back(columns.event_id[i]);
        sample.truth_En.push_back(columns.truth_En[i]);
        sample.truth_omega_n_x.push_back(columns.truth_omega_n_x[i]);
        sample.truth_omega_n_y.push_back(columns.truth_omega_n_y[i]);
        sample.truth_omega_n_z.push_back(columns.truth_omega_n_z[i]);
        sample.truth_Ep.push_back(columns.truth_Ep[i]);
        sample.truth_theta_p.push_back(columns.truth_theta_p[i]);
        sample.truth_p_dir_x.push_back(columns.truth_p_dir_x[i]);
        sample.truth_p_dir_y.push_back(columns.truth_p_dir_y[i]);
        sample.truth_p_dir_z.push_back(columns.truth_p_dir_z[i]);
        sample.det_Evis.push_back(columns.det_Evis[i]);
    }
    return sample;
}

} // namespace riptide::io
