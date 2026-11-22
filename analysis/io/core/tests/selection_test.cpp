#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "riptide/common/event_schema.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"

namespace {

using riptide::io::Cut;
using riptide::io::CutFlow;
using riptide::io::CutInputs;
using riptide::io::IoError;
using riptide::io::Selection;
using riptide::schema::EventClass;

constexpr auto kSingleNp = static_cast<std::int32_t>(EventClass::SingleNp);
constexpr auto kMultipleNp = static_cast<std::int32_t>(EventClass::MultipleNp);

} // namespace

TEST_CASE("selections parse in order and print back", "[selection]") {
    const auto selection = Selection::parse("single_np,contained,evis_min=0.5");
    REQUIRE(selection.cuts().size() == 3);
    CHECK(selection.cuts()[0].cut == Cut::SingleNp);
    CHECK(selection.cuts()[1].cut == Cut::Contained);
    CHECK(selection.cuts()[2].cut == Cut::EvisMin);
    CHECK_THAT(selection.cuts()[2].threshold, Catch::Matchers::WithinULP(0.5, 0));
    CHECK(selection.toString() == "single_np,contained,evis_min=0.5");
    CHECK(Selection::parse(selection.toString()).toString() == selection.toString());
}

TEST_CASE("has_proton parses and reads truth_valid_np", "[selection]") {
    const auto selection = Selection::parse("single_np,has_proton,contained");
    REQUIRE(selection.cuts().size() == 3);
    CHECK(selection.cuts()[1].cut == Cut::HasProton);
    CHECK(selection.toString() == "single_np,has_proton,contained");
    CHECK(riptide::io::cutName(Cut::HasProton) == "has_proton");
    CHECK_FALSE(selection.firstFailingCut({kSingleNp, 0, 2.0, 1}).has_value());
    CHECK(selection.firstFailingCut({kSingleNp, 0, 2.0, 0}) == std::optional<std::size_t>{1});
    // The chain order decides which cut an event fails first.
    CHECK(selection.firstFailingCut({kMultipleNp, 0, 2.0, 0}) == std::optional<std::size_t>{0});
    CHECK(selection.firstFailingCut({kSingleNp, 1, 2.0, 1}) == std::optional<std::size_t>{2});
    CHECK_THROWS_AS(Selection::parse("has_proton=1"), IoError);
    CHECK_THROWS_AS(Selection::parse("has_proton,has_proton"), IoError);
}

TEST_CASE("malformed selections are refused", "[selection]") {
    for (const char* text :
         {"", "single_np,", ",contained", "single_np,single_np", "evis_min",
          "evis_min=", "evis_min=-1", "evis_min=0.5x", "evis_min=nan", "contained=1", "two_np"}) {
        INFO(text);
        CHECK_THROWS_AS(Selection::parse(text), IoError);
    }
}

TEST_CASE("each cut reads its own column", "[selection]") {
    const auto selection = Selection::parse("single_np,contained,evis_min=1");
    CHECK_FALSE(selection.firstFailingCut({kSingleNp, 0, 2.0}).has_value());
    CHECK(selection.firstFailingCut({kMultipleNp, 0, 2.0}) == std::optional<std::size_t>{0});
    CHECK(selection.firstFailingCut({kSingleNp, 1, 2.0}) == std::optional<std::size_t>{1});
    // The threshold is strict: det_Evis must exceed it.
    CHECK(selection.firstFailingCut({kSingleNp, 0, 1.0}) == std::optional<std::size_t>{2});
}

TEST_CASE("the cut flow is monotone and counts the right events", "[selection]") {
    const auto selection = Selection::parse("single_np,contained,evis_min=1");
    const std::vector<CutInputs> events{
        {kSingleNp, 0, 2.0},   // passes all
        {kSingleNp, 0, 0.5},   // fails evis_min
        {kSingleNp, 1, 2.0},   // fails contained
        {kMultipleNp, 0, 2.0}, // fails single_np
        {kSingleNp, 0, 3.0},   // passes all
        {kMultipleNp, 1, 0.0}, // fails single_np first
    };
    CutFlow flow(selection.cuts().size());
    for (const auto& event : events) {
        flow.record(selection.firstFailingCut(event));
    }
    CHECK(flow.total() == 6);
    CHECK(flow.passed(0) == 4);
    CHECK(flow.passed(1) == 3);
    CHECK(flow.passed(2) == 2);
    CHECK(flow.selected() == 2);
    CHECK_THROWS(flow.record(std::optional<std::size_t>{4}));
}
