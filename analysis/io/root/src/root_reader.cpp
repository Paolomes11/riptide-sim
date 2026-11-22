#include "riptide/io/root_reader.hpp"

#include <ROOT/RDataFrame.hxx>

#include <RtypesCore.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "riptide/common/event_schema.hpp"

namespace riptide::io {

namespace {

// The schema's Int32 and Double columns are written by Geant4 as Int_t and
// Double_t leaves; Take<> must use exactly these types.
static_assert(std::is_same_v<std::int32_t, Int_t>);
static_assert(std::is_same_v<double, Double_t>);
constexpr std::string_view kInt32TypeName = "Int_t";
constexpr std::string_view kDoubleTypeName = "Double_t";

template <typename Frame>
void requireColumn(Frame& frame, const std::vector<std::string>& available, std::string_view name,
                   std::string_view type) {
    const std::string column{name};
    if (std::find(available.begin(), available.end(), column) == available.end()) {
        throw IoError("event files lack column " + column);
    }
    const auto actual = frame.GetColumnType(column);
    if (actual != type) {
        throw IoError("column " + column + " has type " + actual + ", expected " +
                      std::string{type});
    }
}

} // namespace

EventColumns readEventColumns(const RunInfo& run, std::uint64_t max_events) {
    std::vector<std::string> files;
    files.reserve(run.event_files.size());
    for (const auto& path : run.event_files) {
        files.push_back(path.string());
    }
    // No implicit multithreading: the entry order of Take<> is then the file
    // order, and the sort below makes the result independent of it anyway.
    ROOT::RDataFrame frame(std::string{schema::kEventTreeName}, files);

    EventColumns columns;
    const auto available = frame.GetColumnNames();
    columns.forEachInt32Column([&](std::string_view name, const std::vector<std::int32_t>&) {
        requireColumn(frame, available, name, kInt32TypeName);
    });
    columns.forEachDoubleColumn([&](std::string_view name, const std::vector<double>&) {
        requireColumn(frame, available, name, kDoubleTypeName);
    });

    const auto expected = expectedEventCount(run.events, max_events);
    // Negative ids are kept so that sortByEventId reports them.
    auto bounded = frame.Filter(
        [expected](Int_t id) { return id < 0 || static_cast<std::uint64_t>(id) < expected; },
        {"event_id"});

    // Book every column first: RDataFrame then fills all of them in one pass.
    std::vector<ROOT::RDF::RResultPtr<std::vector<Int_t>>> int_results;
    std::vector<ROOT::RDF::RResultPtr<std::vector<Double_t>>> double_results;
    columns.forEachInt32Column([&](std::string_view name, const std::vector<std::int32_t>&) {
        int_results.push_back(bounded.Take<Int_t>(std::string{name}));
    });
    columns.forEachDoubleColumn([&](std::string_view name, const std::vector<double>&) {
        double_results.push_back(bounded.Take<Double_t>(std::string{name}));
    });

    std::size_t next_int = 0;
    columns.forEachInt32Column([&](std::string_view, std::vector<std::int32_t>& column) {
        column = std::move(*int_results[next_int++]);
    });
    std::size_t next_double = 0;
    columns.forEachDoubleColumn([&](std::string_view, std::vector<double>& column) {
        column = std::move(*double_results[next_double++]);
    });

    sortByEventId(columns, expected);
    return columns;
}

} // namespace riptide::io
