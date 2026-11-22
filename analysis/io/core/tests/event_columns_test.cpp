#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstdint>
#include <vector>

#include "riptide/common/event_schema.hpp"
#include "riptide/io/event_columns.hpp"
#include "riptide/io/run_info.hpp"

namespace {

using riptide::io::CutFlow;
using riptide::io::EventColumns;
using riptide::io::IoError;
using riptide::io::Selection;
using riptide::schema::EventClass;

constexpr auto kSingleNp = static_cast<std::int32_t>(EventClass::SingleNp);
constexpr auto kMultipleNp = static_cast<std::int32_t>(EventClass::MultipleNp);

// Columns in the order a two-thread run may deliver them; every double column
// carries a value derived from the id so the permutation can be checked.
EventColumns columnsWithIds(const std::vector<std::int32_t>& ids) {
    EventColumns columns;
    for (const auto id : ids) {
        columns.event_id.push_back(id);
        columns.truth_event_class.push_back(id % 2 == 0 ? kSingleNp : kMultipleNp);
        columns.truth_p_exits.push_back(id == 4 ? 1 : 0);
        columns.forEachDoubleColumn(
            [id](std::string_view, std::vector<double>& column) { column.push_back(10.0 + id); });
    }
    return columns;
}

} // namespace

TEST_CASE("columns are sorted by event_id", "[event_columns]") {
    auto columns = columnsWithIds({3, 0, 4, 1, 2});
    riptide::io::sortByEventId(columns, 5);
    CHECK(columns.event_id == std::vector<std::int32_t>{0, 1, 2, 3, 4});
    CHECK(columns.truth_p_exits == std::vector<std::int32_t>{0, 0, 0, 0, 1});
    columns.forEachDoubleColumn([](std::string_view, const std::vector<double>& column) {
        for (std::size_t i = 0; i < column.size(); ++i) {
            CHECK_THAT(column[i], Catch::Matchers::WithinULP(10.0 + static_cast<double>(i), 0));
        }
    });
}

TEST_CASE("missing repeated or extra events are refused", "[event_columns]") {
    auto repeated = columnsWithIds({0, 1, 1});
    CHECK_THROWS_AS(riptide::io::sortByEventId(repeated, 3), IoError);
    auto gap = columnsWithIds({0, 1, 3});
    CHECK_THROWS_AS(riptide::io::sortByEventId(gap, 3), IoError);
    auto negative = columnsWithIds({-1, 0, 1});
    CHECK_THROWS_AS(riptide::io::sortByEventId(negative, 3), IoError);
    auto short_run = columnsWithIds({0, 1});
    CHECK_THROWS_AS(riptide::io::sortByEventId(short_run, 3), IoError);
    auto ragged = columnsWithIds({0, 1, 2});
    ragged.det_Evis.pop_back();
    CHECK_THROWS_AS(riptide::io::sortByEventId(ragged, 3), IoError);
}

TEST_CASE("max_events bounds the expected count", "[event_columns]") {
    CHECK(riptide::io::expectedEventCount(1000, 0) == 1000);
    CHECK(riptide::io::expectedEventCount(1000, 100) == 100);
    CHECK(riptide::io::expectedEventCount(1000, 5000) == 1000);
}

TEST_CASE("selectEvents keeps the selected truth and counts every event", "[event_columns]") {
    auto columns = columnsWithIds({0, 1, 2, 3, 4, 5});
    riptide::io::sortByEventId(columns, 6);
    const auto selection = Selection::parse("single_np,contained");
    CutFlow flow(selection.cuts().size());
    const auto sample = riptide::io::selectEvents(columns, selection, flow);
    // Even ids are single n-p, id 4 exits: 0 and 2 survive.
    CHECK(sample.event_id == std::vector<std::int32_t>{0, 2});
    REQUIRE(sample.truth_Ep.size() == 2);
    CHECK_THAT(sample.truth_Ep[1], Catch::Matchers::WithinULP(12.0, 0));
    CHECK(sample.truth_Ep_phys.empty());
    CHECK(sample.reco_Ep.empty());
    CHECK(flow.total() == 6);
    CHECK(flow.passed(0) == 3);
    CHECK(flow.passed(1) == 2);
    CutFlow wrong(1);
    CHECK_THROWS_AS(riptide::io::selectEvents(columns, selection, wrong), IoError);
}
