#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sstream>
#include <string>

#include "riptide/sim/hp_recoil.hpp"

namespace {

using Catch::Matchers::WithinRel;
using namespace riptide::sim::checks;

// Target mass of H-1 in the G4NDL 4.7.1 elastic file, in neutron masses.
constexpr double kAwrH1 = 0.999167;

// The parsed value is the text itself, up to rounding.
constexpr double kParseRel = 1.0e-14;

} // namespace

TEST_CASE("parseHpElasticHeader reads the G4NDL first line", "[hp_recoil]") {
    std::istringstream in(" \t1 \t0.999167 \t2\n \t96\n");
    const HpElasticHeader h = parseHpElasticHeader(in);
    CHECK(h.rep_flag == 1);
    CHECK_THAT(h.awr, WithinRel(kAwrH1, kParseRel));
    CHECK(h.frame_flag == kHpFrameCm);
}

TEST_CASE("parseHpElasticHeader refuses what the model does not cover", "[hp_recoil]") {
    for (const std::string text : {"1 0.999167 1", "1 -0.5 2", "0 0.999167 2", "1 0.999167", ""}) {
        std::istringstream in(text);
        CHECK_THROWS_AS(parseHpElasticHeader(in), HpHeaderError);
    }
}
