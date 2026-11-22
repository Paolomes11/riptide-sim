#pragma once

// The columns the analysis reads from the event files, as they come from the
// reader, before sorting and selection. Keeping this step free of ROOT lets
// the ordering, the event-id checks and the selection be tested in CI; the
// ROOT reader only fills the vectors by name.

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "riptide/common/event_sample.hpp"
#include "riptide/io/selection.hpp"

namespace riptide::io {

struct EventColumns {
    std::vector<std::int32_t> event_id;
    std::vector<std::int32_t> truth_event_class;
    std::vector<std::int32_t> truth_p_exits;

    std::vector<double> truth_En;
    std::vector<double> truth_omega_n_x;
    std::vector<double> truth_omega_n_y;
    std::vector<double> truth_omega_n_z;
    std::vector<double> truth_Ep;
    std::vector<double> truth_theta_p;
    std::vector<double> truth_p_dir_x;
    std::vector<double> truth_p_dir_y;
    std::vector<double> truth_p_dir_z;
    std::vector<double> det_Evis;

    [[nodiscard]] std::size_t size() const noexcept {
        return event_id.size();
    }

    // visit(name, column) on every column, with its name in the event schema:
    // the reader binds columns through this list, so a column added here is
    // read without touching the reader.
    template <typename F> void forEachInt32Column(F&& visit) {
        visit(std::string_view{"event_id"}, event_id);
        visit(std::string_view{"truth_event_class"}, truth_event_class);
        visit(std::string_view{"truth_p_exits"}, truth_p_exits);
    }
    template <typename F> void forEachDoubleColumn(F&& visit) {
        visit(std::string_view{"truth_En"}, truth_En);
        visit(std::string_view{"truth_omega_n_x"}, truth_omega_n_x);
        visit(std::string_view{"truth_omega_n_y"}, truth_omega_n_y);
        visit(std::string_view{"truth_omega_n_z"}, truth_omega_n_z);
        visit(std::string_view{"truth_Ep"}, truth_Ep);
        visit(std::string_view{"truth_theta_p"}, truth_theta_p);
        visit(std::string_view{"truth_p_dir_x"}, truth_p_dir_x);
        visit(std::string_view{"truth_p_dir_y"}, truth_p_dir_y);
        visit(std::string_view{"truth_p_dir_z"}, truth_p_dir_z);
        visit(std::string_view{"det_Evis"}, det_Evis);
    }
};

// Events whose event_id is below this bound are read; with max_events = 0
// (all events) the bound is the run's event count.
[[nodiscard]] std::uint64_t expectedEventCount(std::uint64_t events_in_run,
                                               std::uint64_t max_events) noexcept;

// Sorts every column by event_id, so the sample does not depend on how the
// events were split among worker threads and files. Throws IoError unless the
// columns have equal lengths and the ids are exactly 0 .. expected - 1, each
// once: a missing or repeated event means a damaged or mixed-up run.
void sortByEventId(EventColumns& columns, std::uint64_t expected);

// Applies the selection in event_id order, counting every event in cut_flow,
// and returns the truth columns of the selected events. truth_Ep_phys and
// the reco_ columns are left empty for the conversion and response steps.
[[nodiscard]] EventSample selectEvents(const EventColumns& columns, const Selection& selection,
                                       CutFlow& cut_flow);

} // namespace riptide::io
