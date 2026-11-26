#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/stats/calibration_c.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_c.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_c.hpp"

namespace st = riptide::stats;
namespace geo = riptide::geometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr std::uint64_t kSeed = 20260907;
const st::CaseBResolution kResolution{0.10, 0.08};
const st::EnergyDomain kDomain{10.0, 20.0, st::kEnergyGridPoints};

st::CaseCModel toyModel() {
    return {st::energyGrid(kDomain), st::recoilAngleGrid(st::kTrackAngleGridPoints),
            st::Grid{0.0, kPi, st::kTrackObservedGridPoints}, kResolution};
}

void requireSameHyper(const st::HyperparameterSummary& a, const st::HyperparameterSummary& b) {
    REQUIRE_THAT(a.mean, WithinULP(b.mean, 0));
    REQUIRE_THAT(a.std, WithinULP(b.std, 0));
    REQUIRE_THAT(a.median, WithinULP(b.median, 0));
    for (std::size_t l = 0; l < a.intervals.size(); ++l) {
        REQUIRE_THAT(a.intervals[l].lo, WithinULP(b.intervals[l].lo, 0));
        REQUIRE_THAT(a.intervals[l].hi, WithinULP(b.intervals[l].hi, 0));
    }
}

void requireSameDirection(const geo::Direction& a, const geo::Direction& b) {
    REQUIRE_THAT(a.x, WithinULP(b.x, 0));
    REQUIRE_THAT(a.y, WithinULP(b.y, 0));
    REQUIRE_THAT(a.z, WithinULP(b.z, 0));
}

// A result with the given posterior of mu_E and log sigma_E: symmetric
// intervals mean +- {1, 1.645, 1.96} std, the median at the mean.
st::HyperparameterSummary handHyper(double mean, double std) {
    constexpr std::array<double, 3> kHalfWidths{1.0, 1.645, 1.96};
    st::HyperparameterSummary h{mean, std, mean, {}};
    for (std::size_t l = 0; l < h.intervals.size(); ++l) {
        h.intervals[l] = {mean - (kHalfWidths[l] * std), mean + (kHalfWidths[l] * std)};
    }
    return h;
}

st::CaseCExperimentResult handResult(double mu_true, double mu_hat, double ls_hat) {
    st::CaseCExperimentResult r;
    r.truth = {mu_true, 1.0, geo::Direction{0.0, 0.0, 1.0}};
    r.hyper = {handHyper(mu_hat, 0.5), handHyper(ls_hat, 0.25)};
    r.omega_hat = r.truth.omega_n;
    r.direction = {0.0, 0.01, 0.0, true, {true, true, true}};
    return r;
}

} // namespace

TEST_CASE("a Case C row does not depend on the number of threads", "[case_c_calibration]") {
    const st::CaseCModel model = toyModel();
    const std::vector<geo::Direction> coarse = st::sphereGrid(st::kSphereGridPixels);
    const st::CaseCRowSpec row{20,  4, 3, st::kCaseCChecklistTruth, st::CaseCSpectrum::truncated,
                               true};
    const auto one = st::runCaseCRow(model, kDomain, row, kResolution, kSeed, coarse, 1);
    const auto four = st::runCaseCRow(model, kDomain, row, kResolution, kSeed, coarse, 4);
    REQUIRE(one.size() == 4);
    REQUIRE(four.size() == 4);
    for (std::size_t k = 0; k < one.size(); ++k) {
        REQUIRE_THAT(one[k].truth.mu, WithinULP(four[k].truth.mu, 0));
        REQUIRE_THAT(one[k].truth.sigma, WithinULP(four[k].truth.sigma, 0));
        requireSameHyper(one[k].hyper.mu, four[k].hyper.mu);
        requireSameHyper(one[k].hyper.log_sigma, four[k].hyper.log_sigma);
        REQUIRE(one[k].alternate.has_value());
        REQUIRE(four[k].alternate.has_value());
        requireSameHyper(one[k].alternate->mu, four[k].alternate->mu);
        requireSameHyper(one[k].alternate->log_sigma, four[k].alternate->log_sigma);
        REQUIRE_THAT(one[k].predictive.mean, WithinULP(four[k].predictive.mean, 0));
        REQUIRE_THAT(one[k].predictive.std, WithinULP(four[k].predictive.std, 0));
        requireSameDirection(one[k].omega_0, four[k].omega_0);
        requireSameDirection(one[k].omega_hat, four[k].omega_hat);
        REQUIRE_THAT(one[k].direction.error, WithinULP(four[k].direction.error, 0));
        REQUIRE_THAT(one[k].direction.sigma, WithinULP(four[k].direction.sigma, 0));
    }
    // Sanity at N = 20: mu_E within a few posterior stds of the truth, and
    // the stage-1' direction within half a radian.
    for (const st::CaseCExperimentResult& r : one) {
        CHECK(std::abs(r.hyper.mu.mean - r.truth.mu) < 5.0 * r.hyper.mu.std);
        CHECK(r.direction.error < 0.5);
    }
    const st::CaseCEnsemble e = st::summariseCaseCEnsemble(one);
    CHECK(e.mu.entries == 4);
    CHECK(e.direction.entries == 4);
    CHECK(st::caseCPriorShift(one).mu.max_abs >= 0.0);
}

TEST_CASE("Case C hyperparameter ensembles follow the toy definitions", "[case_c_calibration]") {
    // mu residuals +0.5, -0.5, +1.0 with std 0.5: pulls +1, -1, +2.
    std::vector<st::CaseCExperimentResult> r{handResult(15.0, 15.5, std::log(1.0) + 0.25),
                                             handResult(14.0, 13.5, std::log(1.0) - 0.25),
                                             handResult(16.0, 17.0, std::log(1.0))};
    // A median off the mean for log sigma: the residual uses the median.
    r[2].hyper.log_sigma.median = 0.1;
    const st::CaseCEnsemble e = st::summariseCaseCEnsemble(r);
    CHECK(e.mu.entries == 3);
    CHECK_THAT(e.mu.bias, WithinRel(1.0 / 3.0, 1e-14));
    CHECK_THAT(e.mu.rms, WithinRel(std::sqrt(1.5 / 3.0), 1e-14));
    CHECK_THAT(e.mu.std_mean, WithinRel(0.5, 1e-14));
    CHECK_THAT(e.mu.pull_mean, WithinRel(2.0 / 3.0, 1e-14));
    // Population std (ddof = 0) of the pulls {1, -1, 2}.
    const double m = 2.0 / 3.0;
    const double width = std::sqrt(
        (((1.0 - m) * (1.0 - m)) + ((-1.0 - m) * (-1.0 - m)) + ((2.0 - m) * (2.0 - m))) / 3.0);
    CHECK_THAT(e.mu.pull_width, WithinRel(width, 1e-14));
    CHECK_THAT(e.mu.pull_mean_error,
               WithinRel(width * std::sqrt(3.0 / 2.0) / std::sqrt(3.0), 1e-14));
    // Pulls 1, -1, 2 against half-widths 1, 1.645, 1.96: 2/3, 2/3, 2/3.
    for (const double c : e.mu.coverage) {
        CHECK_THAT(c, WithinRel(2.0 / 3.0, 1e-14));
    }
    // log sigma: residuals +0.25, -0.25, +0.1 (median); pulls from the means.
    CHECK_THAT(e.log_sigma.bias, WithinRel(0.1 / 3.0, 1e-12));
    CHECK_THAT(e.log_sigma.pull_mean, WithinAbs(0.0, 1e-14));
    CHECK(e.direction.entries == 3);
    CHECK_THAT(e.direction.truth_in_cap_fraction, WithinRel(1.0, 1e-14));
}

TEST_CASE("the Case C prior shift compares estimates on the same experiments",
          "[case_c_calibration]") {
    std::vector<st::CaseCExperimentResult> r{handResult(15.0, 15.0, 0.0),
                                             handResult(15.0, 15.0, 0.0)};
    // Alternate moves mu by +0.25 and -0.5 (std 0.5), log sigma median +0.05
    // in both (std 0.25).
    r[0].alternate = st::CaseCHyperSummary{handHyper(15.25, 0.5), handHyper(0.05, 0.25)};
    r[1].alternate = st::CaseCHyperSummary{handHyper(14.5, 0.5), handHyper(0.05, 0.25)};
    const st::CaseCPriorShift s = st::caseCPriorShift(r);
    CHECK_THAT(s.mu.mean, WithinRel(-0.25, 1e-14));
    CHECK_THAT(s.mu.max_abs, WithinRel(1.0, 1e-14));
    CHECK_THAT(s.mu.mean_error, WithinRel(0.75, 1e-14));
    CHECK_THAT(s.log_sigma.mean, WithinRel(0.2, 1e-12));
    CHECK_THAT(s.log_sigma.mean_error, WithinAbs(0.0, 1e-14));
    const st::CaseCEnsemble alt = st::summariseCaseCEnsemble(r, true);
    CHECK_THAT(alt.mu.bias, WithinRel(-0.125, 1e-14));
}

TEST_CASE("Case C ensembles reject what they cannot summarise", "[case_c_calibration]") {
    const std::vector<st::CaseCExperimentResult> one{handResult(15.0, 15.0, 0.0)};
    CHECK_THROWS_AS(st::summariseCaseCEnsemble(one), std::invalid_argument);
    CHECK_THROWS_AS(st::caseCPriorShift(one), std::invalid_argument);
    const std::vector<st::CaseCExperimentResult> two{handResult(15.0, 15.0, 0.0),
                                                     handResult(15.0, 15.0, 0.0)};
    CHECK_THROWS_AS(st::summariseCaseCEnsemble(two, true), std::invalid_argument);
    CHECK_THROWS_AS(st::caseCPriorShift(two), std::invalid_argument);
}
