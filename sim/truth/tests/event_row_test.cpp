#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string_view>
#include <variant>

#include "riptide/common/units.hpp"
#include "riptide/sim/event_row.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using riptide::schema::ColumnType;
using riptide::schema::FirstInteraction;
using riptide::sim::columnIndex;
using riptide::sim::columnSpec;
using riptide::sim::EventMeta;
using riptide::sim::eventRow;
using riptide::sim::EventRow;
using riptide::sim::kEventColumns;
using riptide::sim::TruthRecord;
using riptide::sim::truthRecord;
namespace u = riptide::units;

TruthRecord knownRecord() {
    TruthRecord r;
    r.En = 14.0 * u::MeV;
    r.omega_n = {0.0, 0.0, 1.0};
    r.r0 = {1.0 * u::mm, -2.0 * u::mm, 3.0 * u::mm};
    r.t0 = 0.25 * u::ns;
    r.first_interaction = FirstInteraction::NCElastic;
    r.n_np_scatters = 2;
    r.Ep = 5.5 * u::MeV;
    r.p_dir = {0.6, 0.0, 0.8};
    r.theta_p = 0.6435 * u::rad;
    r.valid_np = true;
    r.n_np_without_proton = 1;
    return r;
}

double asDouble(const EventRow& row, std::string_view name) {
    return std::get<double>(row.at(columnIndex(name)));
}

std::int32_t asInt(const EventRow& row, std::string_view name) {
    return std::get<std::int32_t>(row.at(columnIndex(name)));
}

} // namespace

TEST_CASE("event row: every value has the schema's column type", "[row]") {
    const EventRow row = eventRow(EventMeta{}, TruthRecord{});
    for (std::size_t i = 0; i < kEventColumns; ++i) {
        INFO("column " << columnSpec(i).name);
        const bool is_int = std::holds_alternative<std::int32_t>(row.at(i));
        CHECK(is_int == (columnSpec(i).type == ColumnType::Int32));
    }
}

TEST_CASE("event row: column names are unique and indexed in schema order", "[row]") {
    std::set<std::string_view> names;
    for (std::size_t i = 0; i < kEventColumns; ++i) {
        names.insert(columnSpec(i).name);
        CHECK(columnIndex(columnSpec(i).name) == i);
    }
    CHECK(names.size() == kEventColumns);
    CHECK(columnIndex("event_id") == 0);
    CHECK_THROWS_AS(columnIndex("not_a_column"), std::out_of_range);
}

TEST_CASE("event row: values land in the column of the same name", "[row]") {
    const EventMeta meta{42, 123456, 20260907};
    const EventRow row = eventRow(meta, knownRecord());
    constexpr double kExact = 0.0;

    CHECK(asInt(row, "event_id") == 42);
    CHECK(asInt(row, "rng_seed_0") == 123456);
    CHECK(asInt(row, "rng_seed_1") == 20260907);
    CHECK_THAT(asDouble(row, "truth_En"), WithinAbs(14.0, kExact));
    CHECK_THAT(asDouble(row, "truth_omega_n_z"), WithinAbs(1.0, kExact));
    CHECK_THAT(asDouble(row, "truth_r0_x"), WithinAbs(1.0, kExact));
    CHECK_THAT(asDouble(row, "truth_r0_y"), WithinAbs(-2.0, kExact));
    CHECK_THAT(asDouble(row, "truth_r0_z"), WithinAbs(3.0, kExact));
    CHECK_THAT(asDouble(row, "truth_t0"), WithinAbs(0.25, kExact));
    CHECK(asInt(row, "truth_first_interaction") ==
          static_cast<std::int32_t>(FirstInteraction::NCElastic));
    CHECK(asInt(row, "truth_n_np_scatters") == 2);
    CHECK_THAT(asDouble(row, "truth_Ep"), WithinAbs(5.5, kExact));
    CHECK_THAT(asDouble(row, "truth_p_dir_x"), WithinAbs(0.6, kExact));
    CHECK_THAT(asDouble(row, "truth_p_dir_y"), WithinAbs(0.0, kExact));
    CHECK_THAT(asDouble(row, "truth_p_dir_z"), WithinAbs(0.8, kExact));
    CHECK_THAT(asDouble(row, "truth_theta_p"), WithinAbs(0.6435, kExact));
    CHECK(asInt(row, "truth_valid_np") == 1);
}

TEST_CASE("event row: valid_np false is written as 0", "[row]") {
    TruthRecord r = knownRecord();
    r.valid_np = false;
    CHECK(asInt(eventRow(EventMeta{}, r), "truth_valid_np") == 0);
}

TEST_CASE("truth record: reading a row back gives the written record", "[row]") {
    const TruthRecord in = knownRecord();
    const TruthRecord out = truthRecord(eventRow(EventMeta{7, 1, 2}, in));
    constexpr double kExact = 0.0;
    CHECK_THAT(out.En, WithinAbs(in.En, kExact));
    CHECK_THAT(out.omega_n.z, WithinAbs(in.omega_n.z, kExact));
    CHECK_THAT(out.r0.x, WithinAbs(in.r0.x, kExact));
    CHECK_THAT(out.r0.y, WithinAbs(in.r0.y, kExact));
    CHECK_THAT(out.r0.z, WithinAbs(in.r0.z, kExact));
    CHECK_THAT(out.t0, WithinAbs(in.t0, kExact));
    CHECK(out.first_interaction == in.first_interaction);
    CHECK(out.n_np_scatters == in.n_np_scatters);
    CHECK_THAT(out.Ep, WithinAbs(in.Ep, kExact));
    CHECK_THAT(out.p_dir.x, WithinAbs(in.p_dir.x, kExact));
    CHECK_THAT(out.p_dir.z, WithinAbs(in.p_dir.z, kExact));
    CHECK_THAT(out.theta_p, WithinAbs(in.theta_p, kExact));
    CHECK(out.valid_np);
    // Not a schema column, so it cannot survive the round trip.
    CHECK(out.n_np_without_proton == 0);
}

TEST_CASE("truth record: every first-interaction code is read back", "[row]") {
    for (const FirstInteraction f :
         {FirstInteraction::Unknown, FirstInteraction::NpElastic, FirstInteraction::NCElastic,
          FirstInteraction::NCInelastic, FirstInteraction::Other}) {
        TruthRecord r;
        r.first_interaction = f;
        CHECK(truthRecord(eventRow(EventMeta{}, r)).first_interaction == f);
    }
}

TEST_CASE("truth record: malformed rows are rejected", "[row]") {
    EventRow row = eventRow(EventMeta{}, knownRecord());
    SECTION("unknown first-interaction code") {
        row.at(columnIndex("truth_first_interaction")) = std::int32_t{5};
        CHECK_THROWS_AS(truthRecord(row), std::out_of_range);
    }
    SECTION("valid_np neither 0 nor 1") {
        row.at(columnIndex("truth_valid_np")) = std::int32_t{2};
        CHECK_THROWS_AS(truthRecord(row), std::out_of_range);
    }
    SECTION("value of the wrong type") {
        row.at(columnIndex("truth_En")) = std::int32_t{14};
        CHECK_THROWS_AS(truthRecord(row), std::bad_variant_access);
    }
}
