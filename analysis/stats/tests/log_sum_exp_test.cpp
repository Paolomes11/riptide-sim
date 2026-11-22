#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "riptide/stats/log_sum_exp.hpp"

namespace st = riptide::stats;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double kInf = std::numeric_limits<double>::infinity();
} // namespace

TEST_CASE("logSumExp matches the direct sum on moderate values", "[stats][logsumexp]") {
    const std::vector<double> values{-1.0, 0.5, 2.0};
    const double direct = std::log(std::exp(-1.0) + std::exp(0.5) + std::exp(2.0));
    CHECK_THAT(st::logSumExp(values), WithinRel(direct, 1e-15));
}

TEST_CASE("logSumExp is stable on extreme values", "[stats][logsumexp]") {
    // exp(-1e4) underflows to 0: the direct sum would give -inf.
    const std::vector<double> tiny{-1e4, -1e4};
    CHECK_THAT(st::logSumExp(tiny), WithinRel(-1e4 + std::log(2.0), 1e-15));
    // exp(1e3) overflows to +inf.
    const std::vector<double> huge{1e3, 1e3, 1e3};
    CHECK_THAT(st::logSumExp(huge), WithinRel(1e3 + std::log(3.0), 1e-15));
    // A term 800 below the maximum contributes nothing at double precision.
    const std::vector<double> spread{0.0, -800.0};
    CHECK_THAT(st::logSumExp(spread), WithinAbs(0.0, 1e-300));
}

TEST_CASE("logSumExp accepts -inf terms", "[stats][logsumexp]") {
    const std::vector<double> some{-kInf, 1.0, -kInf};
    CHECK_THAT(st::logSumExp(some), WithinRel(1.0, 1e-15));
    const std::vector<double> all{-kInf, -kInf};
    const double result = st::logSumExp(all);
    CHECK((std::isinf(result) && result < 0.0));
    const double empty = st::logSumExp(std::vector<double>{});
    CHECK((std::isinf(empty) && empty < 0.0));
}

TEST_CASE("logSumExp rejects NaN and +inf", "[stats][logsumexp]") {
    CHECK_THROWS_AS(st::logSumExp(std::vector<double>{0.0, std::nan("")}), std::invalid_argument);
    CHECK_THROWS_AS(st::logSumExp(std::vector<double>{0.0, kInf}), std::invalid_argument);
}
