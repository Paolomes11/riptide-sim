#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <string>

#include "riptide/sim/endf_format.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace endf = riptide::sim::endf;
using endf::EndfError;
using endf::InterpolationLaw;

TEST_CASE("Fortran numbers without the exponent letter", "[endf][format]") {
    // Fields copied from the H-1 evaluations.
    CHECK_THAT(endf::parseEndfFloat(" 9.991673-1"), WithinULP(9.991673e-1, 0));
    CHECK_THAT(endf::parseEndfFloat(" 2.000000+7"), WithinULP(2.0e7, 0));
    CHECK_THAT(endf::parseEndfFloat("-1.06436-2 "), WithinULP(-1.06436e-2, 0));
    CHECK_THAT(endf::parseEndfFloat("1.001000+3"), WithinULP(1001.0, 0));
    CHECK_THAT(endf::parseEndfFloat("-1.2+3"), WithinULP(-1200.0, 0));
    CHECK_THAT(endf::parseEndfFloat("+4.5-10"), WithinULP(4.5e-10, 0));
}

TEST_CASE("other accepted number forms", "[endf][format]") {
    CHECK_THAT(endf::parseEndfFloat(" 1.0E+00"), WithinULP(1.0, 0));
    CHECK_THAT(endf::parseEndfFloat("2.5e-3"), WithinULP(2.5e-3, 0));
    CHECK_THAT(endf::parseEndfFloat("         42"), WithinULP(42.0, 0));
    CHECK_THAT(endf::parseEndfFloat("-0.5"), WithinULP(-0.5, 0));
    CHECK_THAT(endf::parseEndfFloat("           "), WithinULP(0.0, 0));
    CHECK_THAT(endf::parseEndfFloat(""), WithinULP(0.0, 0));
}

TEST_CASE("malformed numbers are refused", "[endf][format]") {
    const std::string field =
        GENERATE("abc", "1.0-", "1.0 -3", "1.0--3", "1.0+3x", "--1", "1e999", "nan", "1.0,5");
    INFO("field '" << field << "'");
    CHECK_THROWS_AS(endf::parseEndfFloat(field), EndfError);
}

TEST_CASE("integer fields", "[endf][format]") {
    CHECK(endf::parseEndfInt("          2") == 2);
    CHECK(endf::parseEndfInt("        -96") == -96);
    CHECK(endf::parseEndfInt("           ") == 0);
    CHECK_THROWS_AS(endf::parseEndfInt("        2.0"), EndfError);
    CHECK_THROWS_AS(endf::parseEndfInt("       1 2"), EndfError);
}

TEST_CASE("an ENDF line splits into fields and MAT/MF/MT", "[endf][format]") {
    // MF4/MT2 CONT record of the JEFF-3.3 H-1 file.
    const std::string line =
        " 0.000000+0 9.991673-1          0          2          0          0 125 4  2    2";
    const endf::EndfLine parsed = endf::parseEndfLine(line);
    CHECK(parsed.mat == 125);
    CHECK(parsed.mf == 4);
    CHECK(parsed.mt == 2);
    CHECK_THAT(parsed.real(1), WithinULP(9.991673e-1, 0));
    CHECK(parsed.integer(3) == 2);
    // The sequence number is optional, and a trailing CR is tolerated.
    const endf::EndfLine short_line = endf::parseEndfLine(line.substr(0, 75) + "\r");
    CHECK(short_line.mt == 2);
    CHECK_THROWS_AS(endf::parseEndfLine(line.substr(0, 74)), EndfError);
}

TEST_CASE("interpolation law codes", "[endf][interpolation]") {
    CHECK(endf::interpolationLaw(2) == InterpolationLaw::LinLin);
    CHECK(endf::interpolationLaw(5) == InterpolationLaw::LogLog);
    CHECK_THROWS_AS(endf::interpolationLaw(0), EndfError);
    CHECK_THROWS_AS(endf::interpolationLaw(6), EndfError);
}

TEST_CASE("each law is exact on its own functional form", "[endf][interpolation]") {
    const double x1 = 2.0;
    const double x2 = 8.0;
    const double x = GENERATE(2.0, 3.1, 5.0, 7.99, 8.0);
    INFO("x = " << x);
    // 1: constant at the lower value.
    CHECK_THAT(endf::interpolate(InterpolationLaw::Histogram, x1, 3.0, x2, 9.0, x),
               WithinULP(3.0, 0));
    // 2: straight line y = 1.5 x - 2.
    const auto line = [](double v) { return 1.5 * v - 2.0; };
    CHECK_THAT(endf::interpolate(InterpolationLaw::LinLin, x1, line(x1), x2, line(x2), x),
               WithinRel(line(x), 1e-14));
    // 3: y = 0.7 ln(x) + 1.
    const auto lin_log = [](double v) { return 0.7 * std::log(v) + 1.0; };
    CHECK_THAT(endf::interpolate(InterpolationLaw::LinLog, x1, lin_log(x1), x2, lin_log(x2), x),
               WithinRel(lin_log(x), 1e-14));
    // 4: y = 3 exp(-0.4 x).
    const auto log_lin = [](double v) { return 3.0 * std::exp(-0.4 * v); };
    CHECK_THAT(endf::interpolate(InterpolationLaw::LogLin, x1, log_lin(x1), x2, log_lin(x2), x),
               WithinRel(log_lin(x), 1e-13));
    // 5: power law y = 5 x^-0.5 (the 1/v shape of capture cross sections).
    const auto log_log = [](double v) { return 5.0 / std::sqrt(v); };
    CHECK_THAT(endf::interpolate(InterpolationLaw::LogLog, x1, log_log(x1), x2, log_log(x2), x),
               WithinRel(log_log(x), 1e-13));
}

TEST_CASE("logarithmic laws refuse non-positive values", "[endf][interpolation]") {
    CHECK_THROWS_AS(endf::interpolate(InterpolationLaw::LinLog, 0.0, 1.0, 1.0, 2.0, 0.5),
                    EndfError);
    CHECK_THROWS_AS(endf::interpolate(InterpolationLaw::LogLin, 0.0, 0.0, 1.0, 2.0, 0.5),
                    EndfError);
    CHECK_THROWS_AS(endf::interpolate(InterpolationLaw::LogLog, 1.0, 1.0, 2.0, -1.0, 1.5),
                    EndfError);
}

namespace {

// Two ranges: log-log on points 1-3, lin-lin on points 3-5 (NBT = 3, 5), with
// the shape of Figure 1 of ENDF-102 in miniature.
endf::Tab1 twoRangeTable() {
    endf::Tab1 table;
    table.boundaries = {3, 5};
    table.laws = {InterpolationLaw::LogLog, InterpolationLaw::LinLin};
    table.x = {1.0, 4.0, 16.0, 20.0, 30.0};
    table.y = {8.0, 4.0, 2.0, 6.0, 1.0};
    return table;
}

} // namespace

TEST_CASE("TAB1 evaluation picks the range of each interval", "[endf][interpolation]") {
    const endf::Tab1 table = twoRangeTable();
    REQUIRE_NOTHROW(endf::validate(table, "test"));
    // Tabulated points are returned exactly.
    for (std::size_t i = 0; i < table.x.size(); ++i) {
        CHECK_THAT(endf::evaluate(table, table.x[i]), WithinULP(table.y[i], 0));
    }
    // Range 1 (log-log, y = 8 x^-1/2): x = 9 lies in the interval (4, 16).
    CHECK_THAT(endf::evaluate(table, 9.0), WithinRel(8.0 / 3.0, 1e-14));
    // The interval (16, 20) ends at point 4 <= NBT(2): lin-lin, not log-log.
    CHECK_THAT(endf::evaluate(table, 18.0), WithinRel(4.0, 1e-14));
    CHECK_THAT(endf::evaluate(table, 25.0), WithinRel(3.5, 1e-14));
}

TEST_CASE("TAB1 discontinuities and range limits", "[endf][interpolation]") {
    endf::Tab1 table;
    table.boundaries = {4};
    table.laws = {InterpolationLaw::LinLin};
    table.x = {0.0, 1.0, 1.0, 2.0};
    table.y = {0.0, 1.0, 5.0, 7.0};
    REQUIRE_NOTHROW(endf::validate(table, "test"));
    CHECK_THAT(endf::evaluate(table, 0.5), WithinAbs(0.5, 1e-15));
    // At the repeated x the value above the step is used.
    CHECK_THAT(endf::evaluate(table, 1.0), WithinULP(5.0, 0));
    CHECK_THAT(endf::evaluate(table, 1.5), WithinRel(6.0, 1e-15));
    CHECK_THROWS_AS(endf::evaluate(table, -1e-9), EndfError);
    CHECK_THROWS_AS(endf::evaluate(table, 2.0 + 1e-9), EndfError);
    CHECK_THROWS_AS(endf::evaluate(table, std::nan("")), EndfError);
}

TEST_CASE("malformed TAB1 tables are refused", "[endf][interpolation]") {
    endf::Tab1 table = twoRangeTable();
    SECTION("boundary not ending at the last point") {
        table.boundaries = {3, 4};
    }
    SECTION("boundaries not increasing") {
        table.boundaries = {3, 3};
    }
    SECTION("law count differs from range count") {
        table.laws.pop_back();
    }
    SECTION("x decreasing") {
        table.x[3] = 10.0;
    }
    SECTION("x repeated three times") {
        table.x = {1.0, 4.0, 4.0, 4.0, 30.0};
    }
    SECTION("y missing") {
        table.y.pop_back();
    }
    CHECK_THROWS_AS(endf::validate(table, "test"), EndfError);
}
