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

TEST_CASE("first-interaction codes are stable", "[schema]") {
    using riptide::schema::FirstInteraction;
    // Stored in files: changing a value would silently relabel old data.
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::Unknown) == 0);
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::NpElastic) == 1);
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::NCElastic) == 2);
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::NCInelastic) == 3);
    STATIC_REQUIRE(static_cast<int>(FirstInteraction::Other) == 4);
}
