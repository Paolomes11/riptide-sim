#pragma once

// The only ROOT-dependent part of the analysis input: it fills EventColumns
// from the event files of a run through RDataFrame. Everything after this
// (ordering, checks, selection, response) is plain C++ in io_core.

#include <cstdint>

#include "riptide/io/event_columns.hpp"
#include "riptide/io/run_info.hpp"

namespace riptide::io {

// Reads every column of EventColumns for the events of `run` with
// event_id < expectedEventCount(run.events, max_events), sorted by event_id.
// Throws IoError if a column is missing or does not have the type of the event
// schema (a column is never replaced by zeros: a silent zero would bias the
// analysis), or if the ids are not exactly 0 .. expected - 1.
[[nodiscard]] EventColumns readEventColumns(const RunInfo& run, std::uint64_t max_events);

} // namespace riptide::io
