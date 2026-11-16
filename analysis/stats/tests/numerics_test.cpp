#include <catch2/catch_test_macros.hpp>

#include "riptide/stats/numerics.hpp"

TEST_CASE("the stats library keeps IEEE NaN/inf semantics", "[stats][numerics]") {
    REQUIRE(riptide::stats::ieeeSemanticsPreserved());
}
