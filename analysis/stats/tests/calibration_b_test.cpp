#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/calibration_b.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"
#include "riptide/stats/split.hpp"

namespace {

namespace st = riptide::stats;
namespace geo = riptide::geometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;
using riptide::units::pi;

constexpr std::uint64_t kSeed = 20260907;
constexpr st::CaseBResolution kToyResolution{0.10, 0.08}; // MeV, rad
constexpr double kWideEnMax = 30.0;                       // MeV

st::CaseBModel modelOn(const st::EnergyDomain& domain) {
    const st::Grid En = st::energyGrid(domain);
    return {En, st::recoilAngleGrid(st::kTrackAngleGridPoints),
            st::Grid{0.0, pi, st::kTrackObservedGridPoints}, st::flatEnergyLogPrior(En),
            kToyResolution};
}

// Bit for bit: every number of one result equals the other's (0 ULP).
void requireSameResult(const st::CaseBExperimentResult& a, const st::CaseBExperimentResult& b) {
    for (const auto& [u, v] :
         {std::pair{a.omega_true, b.omega_true}, std::pair{a.omega_hat, b.omega_hat}}) {
        REQUIRE_THAT(u.x, WithinULP(v.x, 0));
        REQUIRE_THAT(u.y, WithinULP(v.y, 0));
        REQUIRE_THAT(u.z, WithinULP(v.z, 0));
    }
    REQUIRE_THAT(a.summary.error, WithinULP(b.summary.error, 0));
    REQUIRE_THAT(a.summary.sigma, WithinULP(b.summary.sigma, 0));
    REQUIRE_THAT(a.summary.pull, WithinULP(b.summary.pull, 0));
    REQUIRE(a.summary.truth_in_cap == b.summary.truth_in_cap);
    REQUIRE(a.summary.covered == b.summary.covered);
}

st::CaseBExperimentResult handResult(double error, double sigma, double pull, bool in_cap,
                                     bool covered_68) {
    st::CaseBExperimentResult r;
    r.summary.error = error;
    r.summary.sigma = sigma;
    r.summary.pull = pull;
    r.summary.truth_in_cap = in_cap;
    r.summary.covered = {covered_68, true, true};
    return r;
}

} // namespace

TEST_CASE("parallelForEach runs every index once for any number of threads",
          "[stats][calibration_b]") {
    constexpr std::size_t kCount = 7;
    for (const unsigned threads : {1U, 2U, 3U, 4U, 8U}) {
        std::vector<std::atomic<int>> hits(kCount);
        st::parallelForEach(kCount, threads, [&](std::size_t k) { ++hits[k]; });
        for (const auto& h : hits) {
            REQUIRE(h.load() == 1);
        }
    }
    std::atomic<int> calls{0};
    st::parallelForEach(0, 4, [&](std::size_t) { ++calls; });
    REQUIRE(calls.load() == 0);
    REQUIRE_THROWS_AS(st::parallelForEach(3, 0, [](std::size_t) {}), std::invalid_argument);
}

TEST_CASE("parallelForEach rethrows the failure of the lowest index", "[stats][calibration_b]") {
    for (const unsigned threads : {1U, 2U, 4U}) {
        const auto task = [](std::size_t k) {
            if (k == 3 || k == 5) {
                throw std::runtime_error(std::to_string(k));
            }
        };
        REQUIRE_THROWS_WITH(st::parallelForEach(8, threads, task), "3");
    }
}

TEST_CASE("Case B checklist rows give the same numbers on 1 and 4 threads",
          "[stats][calibration_b]") {
    const st::EnergyDomain domain{};
    const st::CaseBModel model = modelOn(domain);
    const auto coarse_grid = st::sphereGrid(st::kSphereGridPixels);
    const st::CaseBRowSpec row{10, 6, 0};

    const auto one = st::runCaseBRow(model, domain, row, kToyResolution, kSeed, coarse_grid, 1);
    const auto four = st::runCaseBRow(model, domain, row, kToyResolution, kSeed, coarse_grid, 4);
    REQUIRE(one.size() == row.experiments);
    REQUIRE(four.size() == row.experiments);
    for (std::size_t k = 0; k < one.size(); ++k) {
        requireSameResult(one[k], four[k]);
    }
    // Experiment k is the pseudo-experiment of (section, k), whatever the row.
    const auto e2 =
        st::caseBPseudoExperiment(domain, row.events, kToyResolution, kSeed, row.section, 2);
    REQUIRE_THAT(one[2].omega_true.x, WithinULP(e2.truth_omega_n.x, 0));

    SECTION("and so do the prior rows, on the very experiments of the row") {
        const st::CaseBModel wide = modelOn(st::sameStepDomain(domain, kWideEnMax));
        const st::CaseBPriorSpec prior{1.0, 0.5 * pi, 1};
        const st::CaseBRowSpec short_row{10, 4, 0};
        const auto p1 = st::runCaseBPriorRow(model, wide, domain, short_row, prior, kToyResolution,
                                             kSeed, coarse_grid, 1);
        const auto p4 = st::runCaseBPriorRow(model, wide, domain, short_row, prior, kToyResolution,
                                             kSeed, coarse_grid, 4);
        for (std::size_t k = 0; k < short_row.experiments; ++k) {
            requireSameResult(p1.vmf[k], p4.vmf[k]);
            requireSameResult(p1.energy[k], p4.energy[k]);
            REQUIRE_THAT(p1.vmf[k].omega_true.x, WithinULP(one[k].omega_true.x, 0));
            REQUIRE_THAT(p1.energy[k].omega_true.z, WithinULP(one[k].omega_true.z, 0));
        }
        REQUIRE_THROWS_AS(st::runCaseBPriorRow(model, wide, domain, short_row, {-1.0, 0.0, 1},
                                               kToyResolution, kSeed, coarse_grid, 1),
                          std::invalid_argument);
        REQUIRE_THROWS_AS(st::runCaseBPriorRow(model, wide, domain, short_row, {1.0, 4.0, 1},
                                               kToyResolution, kSeed, coarse_grid, 1),
                          std::invalid_argument);
    }
}

TEST_CASE("Case B pool rows analyse the split experiments, the same on 1 and 4 threads",
          "[stats][calibration_b]") {
    const st::EnergyDomain domain{};
    const st::CaseBModel model = modelOn(domain);
    const auto coarse_grid = st::sphereGrid(st::kSphereGridPixels);
    const geo::Direction omega = geo::normalized({3.0, 2.0, 1.0});
    constexpr std::size_t kPool = 100;
    constexpr std::size_t kEvents = 10;
    constexpr std::size_t kExperiments = 6;
    const auto pool = st::caseBPseudoSample(domain, omega, kPool, kToyResolution, kSeed);
    auto engine = riptide::random::perEventEngine(kSeed, 0);
    const auto split = st::disjointExperiments(kPool, kEvents, kExperiments, engine);

    const auto one = st::runCaseBPoolRow(model, pool, omega, split, coarse_grid, 1);
    const auto four = st::runCaseBPoolRow(model, pool, omega, split, coarse_grid, 4);
    REQUIRE(one.size() == kExperiments);
    REQUIRE(four.size() == kExperiments);
    for (std::size_t k = 0; k < kExperiments; ++k) {
        requireSameResult(one[k], four[k]);
        // Experiment k is the analysis of exactly the events at split[k].
        const st::CaseBPseudoExperiment by_hand{omega, st::caseBPseudoSampleSubset(pool, split[k])};
        requireSameResult(one[k], st::analyseCaseBExperiment(model, by_hand, coarse_grid,
                                                             st::uniformDirectionLogPrior));
        REQUIRE_THAT(one[k].omega_true.x, WithinULP(omega.x, 0));
    }
    REQUIRE(st::runCaseBPoolRow(model, pool, omega, {}, coarse_grid, 4).empty());
    REQUIRE_THROWS_AS(st::runCaseBPoolRow(model, pool, omega, {{kPool}}, coarse_grid, 1),
                      std::invalid_argument);
}

TEST_CASE("summariseCaseBEnsemble on known entries", "[stats][calibration_b]") {
    // Four entries; sqrt(2) pull = sqrt(2) {1, 2.5, 1.5, 2}, so one tail entry
    // (sqrt(2) 2.5 ~ 3.54 > 3) and median sqrt(2) 1.75. Tolerance 1e-14
    // relative: a few roundings on O(1) sums.
    const std::vector<st::CaseBExperimentResult> r{
        handResult(0.1, 0.1, 1.0, true, true), handResult(0.2, 0.1, 2.5, true, false),
        handResult(0.3, 0.2, 1.5, true, true), handResult(0.4, 0.2, 2.0, false, false)};
    const st::CaseBEnsemble e = st::summariseCaseBEnsemble(r);
    constexpr double kRel = 1e-14;
    REQUIRE(e.entries == 4);
    REQUIRE_THAT(e.error_mean, WithinRel(0.25, kRel));
    REQUIRE_THAT(e.error_mean_error, WithinRel(std::sqrt(0.05 / 12.0), kRel));
    REQUIRE_THAT(e.error_rms, WithinRel(std::sqrt(0.075), kRel));
    REQUIRE_THAT(e.sigma_mean, WithinRel(0.15, kRel));
    REQUIRE_THAT(e.sigma_mean_error, WithinRel(std::sqrt(0.01 / 3.0 / 4.0), kRel));
    REQUIRE_THAT(e.pull_rms, WithinRel(std::sqrt(13.5 / 4.0), kRel));
    REQUIRE_THAT(e.scaled_pull_median, WithinRel(std::sqrt(2.0) * 1.75, kRel));
    REQUIRE_THAT(e.scaled_pull_tail_fraction, WithinRel(0.25, kRel));
    REQUIRE_THAT(e.truth_in_cap_fraction, WithinRel(0.75, kRel));
    REQUIRE_THAT(e.coverage[0], WithinRel(0.5, kRel));
    REQUIRE_THAT(e.coverage[2], WithinRel(1.0, kRel));

    REQUIRE_THROWS_AS(st::summariseCaseBEnsemble({r[0]}), std::invalid_argument);
}

TEST_CASE("logLogSlope recovers a power law", "[stats][calibration_b]") {
    // rms = 0.25 N^-1/2 exactly: slope -0.5 up to the roundings of log, 1e-12.
    const std::vector<double> n{30.0, 100.0, 300.0, 1000.0};
    std::vector<double> rms;
    for (const double x : n) {
        rms.push_back(0.25 / std::sqrt(x));
    }
    REQUIRE_THAT(st::logLogSlope(n, rms), WithinAbs(-0.5, 1e-12));
    REQUIRE_THROWS_AS(st::logLogSlope({10.0}, {0.1}), std::invalid_argument);
    REQUIRE_THROWS_AS(st::logLogSlope({10.0, 20.0}, {0.1}), std::invalid_argument);
    REQUIRE_THROWS_AS(st::logLogSlope({10.0, 20.0}, {0.1, 0.0}), std::invalid_argument);
    REQUIRE_THROWS_AS(st::logLogSlope({10.0, 10.0}, {0.1, 0.2}), std::invalid_argument);
}

TEST_CASE("tangentCoordinates follow the polar basis of the truth", "[stats][calibration_b]") {
    // Truth at theta0 in the x-z plane; moving theta by a is +a along e_theta,
    // moving phi gives a pure e_phi offset whose length is the angle. 1e-12
    // absolute on angles of 0.01 rad.
    constexpr double kTheta0 = 1.0;
    constexpr double kAngle = 0.01;
    const geo::Direction truth{std::sin(kTheta0), 0.0, std::cos(kTheta0)};
    const geo::Direction down{std::sin(kTheta0 + kAngle), 0.0, std::cos(kTheta0 + kAngle)};
    const auto t = st::tangentCoordinates(down, truth);
    REQUIRE_THAT(t[0], WithinAbs(kAngle, 1e-12));
    REQUIRE_THAT(t[1], WithinAbs(0.0, 1e-12));

    const geo::Direction side{std::sin(kTheta0) * std::cos(kAngle),
                              std::sin(kTheta0) * std::sin(kAngle), std::cos(kTheta0)};
    const auto s = st::tangentCoordinates(side, truth);
    REQUIRE(s[1] > 0.0);
    REQUIRE_THAT(std::hypot(s[0], s[1]), WithinAbs(geo::angleBetween(side, truth), 1e-12));

    const auto zero = st::tangentCoordinates(truth, truth);
    REQUIRE_THAT(zero[0], WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(zero[1], WithinAbs(0.0, 1e-12));
}

TEST_CASE("directionPriorShift measures the move of the estimate", "[stats][calibration_b]") {
    // Reference estimates at the truth, alternatives a = 0.01 rad away along
    // e_theta: every scalar shift is a / sigma, the vector mean (a / sigma, 0)
    // with zero error. 1e-9 relative: the angle of two unit vectors 0.01 rad
    // apart from rounded components.
    constexpr double kAngle = 0.01;
    constexpr double kSigma = 0.05;
    std::vector<st::CaseBExperimentResult> ref;
    std::vector<st::CaseBExperimentResult> alt;
    for (const double theta0 : {0.7, 1.2, 2.0}) {
        st::CaseBExperimentResult r = handResult(0.0, kSigma, 0.0, true, true);
        r.omega_true = {std::sin(theta0), 0.0, std::cos(theta0)};
        r.omega_hat = r.omega_true;
        ref.push_back(r);
        r.omega_hat = {std::sin(theta0 + kAngle), 0.0, std::cos(theta0 + kAngle)};
        alt.push_back(r);
    }
    const st::DirectionPriorShift s = st::directionPriorShift(ref, alt);
    REQUIRE_THAT(s.shift_mean, WithinRel(kAngle / kSigma, 1e-9));
    REQUIRE_THAT(s.shift_max, WithinRel(kAngle / kSigma, 1e-9));
    REQUIRE_THAT(s.shift_error, WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(s.vector_mean[0], WithinRel(kAngle / kSigma, 1e-9));
    REQUIRE_THAT(s.vector_mean[1], WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(s.vector_error[0], WithinAbs(0.0, 1e-9));

    const st::DirectionPriorShift none = st::directionPriorShift(ref, ref);
    REQUIRE_THAT(none.shift_mean, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(none.vector_mean[0], WithinAbs(0.0, 1e-12));

    REQUIRE_THROWS_AS(st::directionPriorShift(ref, {alt[0], alt[1]}), std::invalid_argument);
    auto moved = alt;
    moved[1].omega_true = moved[0].omega_true;
    REQUIRE_THROWS_AS(st::directionPriorShift(ref, moved), std::invalid_argument);
}
