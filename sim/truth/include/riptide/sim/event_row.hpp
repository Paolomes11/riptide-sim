#pragma once

// One ntuple row of an event, in the column order of event_schema.hpp (meta
// columns, truth columns, then the v2, v3 and v4 columns). Geant4-free, so the mapping from the
// in-memory record to the schema is unit-tested on its own; the Geant4 writer only iterates over
// the row.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <variant>

#include "riptide/common/event_schema.hpp"
#include "riptide/sim/deposit_record.hpp"
#include "riptide/sim/truth_record.hpp"

namespace riptide::sim {

// Value of one column; the alternative matches schema::ColumnType (Int32 first).
using ColumnValue = std::variant<std::int32_t, double>;

inline constexpr std::size_t kEventColumns = schema::kV4Columns;

using EventRow = std::array<ColumnValue, kEventColumns>;

// Schema entry of column i of an EventRow.
[[nodiscard]] const schema::BranchSpec& columnSpec(std::size_t i);

// Position of a column in an EventRow; throws std::out_of_range for a name not
// in the schema.
[[nodiscard]] std::size_t columnIndex(std::string_view name);

// Number of leading EventRow columns present in a file of the given schema
// version: later versions only append columns. Throws std::out_of_range for a
// version this build cannot read.
[[nodiscard]] std::size_t columnCount(std::int32_t schema_version);

// The np-vertex and recoil-proton columns are written as 0 unless
// truth.valid_np, whatever the records hold, so that "0 unless valid_np"
// holds by construction for every column of the proton.
[[nodiscard]] EventRow eventRow(const EventMeta& meta, const TruthRecord& truth,
                                const DepositRecord& deposit);

// Inverse of eventRow for the truth part, so a reader rebuilds the record by
// column name instead of by position. n_np_without_proton is not in the
// schema and stays 0. truth_event_class is derived, so it is not read back
// (in a v1 file it is 0 although the event may have interacted). Throws
// std::bad_variant_access for a value of the wrong type and std::out_of_range
// for a first-interaction code outside the schema or a valid_np flag other
// than 0 and 1: such a row cannot come from eventRow. In a file older than v4
// n_start reads as 0.
[[nodiscard]] TruthRecord truthRecord(const EventRow& row);

// Inverse of eventRow for the deposit part. proton.seen is not a column: it is
// rebuilt from truth_valid_np, the condition under which the proton columns
// are filled. In a file older than v3, proton.n_nuclear_recoil reads as 0 because
// the column is missing, not because there were none. Throws std::out_of_range for
// a flag other than 0 and 1.
[[nodiscard]] DepositRecord depositRecord(const EventRow& row);

} // namespace riptide::sim
