#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <vector>

#include "riptide/common/event_sample.hpp"

namespace {

// Number of double columns of EventSample; update with the struct.
constexpr std::size_t kDoubleColumns = 16;

} // namespace

TEST_CASE("an empty sample is consistent", "[event_sample]") {
    const riptide::EventSample sample;
    CHECK(sample.size() == 0);
    CHECK(sample.consistent());
}

TEST_CASE("the visitor reaches every double column once", "[event_sample]") {
    riptide::EventSample sample;
    std::size_t visited = 0;
    sample.forEachDoubleColumn([&](std::vector<double>& column) {
        column.push_back(static_cast<double>(visited));
        ++visited;
    });
    CHECK(visited == kDoubleColumns);
    // Each column got exactly one entry, so none was visited twice.
    sample.event_id.push_back(0);
    CHECK(sample.consistent());
    CHECK(sample.truth_En.size() == 1);
    CHECK(sample.reco_p_dir_z.size() == 1);
    // A column added to the struct but not to the visitor changes the size:
    // EventSample holds only vectors, so it has no padding between them.
    CHECK(sizeof(riptide::EventSample) ==
          sizeof(sample.event_id) + kDoubleColumns * sizeof(std::vector<double>));
}

TEST_CASE("a column of another length makes the sample inconsistent", "[event_sample]") {
    riptide::EventSample sample;
    sample.event_id = {0, 1};
    sample.forEachDoubleColumn([](std::vector<double>& column) { column = {1.0, 2.0}; });
    REQUIRE(sample.consistent());
    sample.reco_theta_p.pop_back();
    CHECK_FALSE(sample.consistent());
}

TEST_CASE("reserve reaches every column", "[event_sample]") {
    constexpr std::size_t kReserve = 64;
    riptide::EventSample sample;
    sample.reserve(kReserve);
    CHECK(sample.event_id.capacity() >= kReserve);
    const riptide::EventSample& view = sample;
    view.forEachDoubleColumn(
        [&](const std::vector<double>& column) { CHECK(column.capacity() >= kReserve); });
}
