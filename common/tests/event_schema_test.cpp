#include <catch2/catch_test_macros.hpp>

#include <set>
#include <string>
#include <string_view>

#include "riptide/common/event_schema.hpp"

namespace {

bool isSnakeCase(std::string_view s) {
    if (s.empty()) {
        return false;
    }
    for (const char c : s) {
        // Upper case is allowed because physics symbols (En, Ep) keep their case.
        const bool ok =
            (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
        if (!ok) {
            return false;
        }
    }
    return true;
}

} // namespace

TEST_CASE("branch names are unique and well formed", "[schema]") {
    using namespace riptide::schema;
    std::set<std::string_view> names;
    for (const auto& b : kMetaBranches) {
        CHECK(isSnakeCase(b.name));
        CHECK(names.insert(b.name).second);
    }
    for (const auto& b : kTruthBranches) {
        CHECK(isSnakeCase(b.name));
        CHECK(names.insert(b.name).second);
    }
    for (const auto& b : kV2Branches) {
        CHECK(isSnakeCase(b.name));
        CHECK(names.insert(b.name).second);
    }
    for (const auto& b : kV3Branches) {
        CHECK(isSnakeCase(b.name));
        CHECK(names.insert(b.name).second);
    }
}

TEST_CASE("truth branches carry the truth_ prefix and a unit", "[schema]") {
    using namespace riptide::schema;
    for (const auto& b : kTruthBranches) {
        INFO(std::string{b.name});
        CHECK(b.name.substr(0, 6) == "truth_");
        CHECK_FALSE(b.unit.empty());
        CHECK_FALSE(b.description.empty());
    }
}

TEST_CASE("v2 branches carry a truth_ or det_ prefix and a unit", "[schema]") {
    using namespace riptide::schema;
    for (const auto& b : kV2Branches) {
        INFO(std::string{b.name});
        CHECK((b.name.substr(0, 6) == "truth_" || b.name.substr(0, 4) == "det_"));
        CHECK_FALSE(b.unit.empty());
        CHECK_FALSE(b.description.empty());
    }
}

TEST_CASE("v3 branches carry the truth_ prefix, a unit and an integer type", "[schema]") {
    using namespace riptide::schema;
    STATIC_REQUIRE(kV3Columns == kV2Columns + kV3Branches.size());
    for (const auto& b : kV3Branches) {
        INFO(std::string{b.name});
        CHECK(b.name.substr(0, 6) == "truth_");
        CHECK(b.type == ColumnType::Int32);
        CHECK_FALSE(b.unit.empty());
        CHECK_FALSE(b.description.empty());
    }
}

TEST_CASE("v4 branches are the neutron start point in mm", "[schema]") {
    using namespace riptide::schema;
    STATIC_REQUIRE(kV4Columns == kV3Columns + kV4Branches.size());
    for (const auto& b : kV4Branches) {
        INFO(std::string{b.name});
        CHECK(b.name.substr(0, 14) == "truth_n_start_");
        CHECK(b.type == ColumnType::Double);
        CHECK(b.unit == "mm");
        CHECK_FALSE(b.description.empty());
    }
}

TEST_CASE("first-interaction codes are stable", "[schema]") {
    using riptide::schema::FirstInteraction;
    // Stored in files: changing a value would silently relabel old data.
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::Unknown) == 0);
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::NpElastic) == 1);
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::NCElastic) == 2);
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::NCInelastic) == 3);
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::Other) == 4);
}

TEST_CASE("event-class codes are stable", "[schema]") {
    using riptide::schema::EventClass;
    STATIC_REQUIRE(static_cast<int>(EventClass::NoInteraction) == 0);
    STATIC_REQUIRE(static_cast<int>(EventClass::SingleNp) == 1);
    STATIC_REQUIRE(static_cast<int>(EventClass::MultipleNp) == 2);
    STATIC_REQUIRE(static_cast<int>(EventClass::NCElasticOnly) == 3);
    STATIC_REQUIRE(static_cast<int>(EventClass::NCInelasticNoNp) == 4);
    STATIC_REQUIRE(static_cast<int>(EventClass::Mixed) == 5);
    STATIC_REQUIRE(static_cast<int>(EventClass::Other) == 6);
}

TEST_CASE("schema uses only the column types Geant4 can write", "[schema]") {
    using namespace riptide::schema;
    STATIC_REQUIRE(kSchemaVersion == 4);
    STATIC_REQUIRE(kOldestReadableSchemaVersion == 1);
    STATIC_REQUIRE(kV2Columns == kV1Columns + kV2Branches.size());
    for (const auto& b : kV2Branches) {
        INFO(std::string{b.name});
        CHECK((b.type == ColumnType::Int32 || b.type == ColumnType::Double));
    }
    for (const auto& b : kMetaBranches) {
        INFO(std::string{b.name});
        CHECK(b.type == ColumnType::Int32);
    }
    for (const auto& b : kTruthBranches) {
        INFO(std::string{b.name});
        CHECK((b.type == ColumnType::Int32 || b.type == ColumnType::Double));
    }
}
