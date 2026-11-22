#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <vector>

#include "riptide/stats/calibration.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/combine.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment.hpp"
#include "riptide/stats/validate.hpp"

namespace st = riptide::stats;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace {

constexpr std::uint64_t kSeed = 20260907;
constexpr st::CaseAResolution kToyResolution{0.10, 0.08};
const std::vector<double> kLevels{0.68, 0.90};

// A summary entry with a symmetric interval of half-width `half` around mean.
void push(st::PosteriorSummaries& s, double truth, double mean, double std, double half) {
    s.truth.push_back(truth);
    s.mean.push_back(mean);
    s.std.push_back(std);
    for (std::size_t l = 0; l < s.intervals.levels.size(); ++l) {
        s.intervals.values.push_back({mean - half, mean + half});
    }
}

struct Setup {
    st::EnergyDomain domain;
    st::Grid En = st::energyGrid(domain);
    st::Grid theta = st::recoilAngleGrid();
    st::CaseAModel model{En, theta, st::physicalRecoilAngleLogPrior(theta), kToyResolution};
    std::vector<double> log_prior = st::flatEnergyLogPrior(En);
};

} // namespace

TEST_CASE("summaries need valid credible levels", "[stats][calibration]") {
    CHECK_THROWS_AS(st::makeSummaries({}), std::invalid_argument);
    CHECK_THROWS_AS(st::makeSummaries({0.68, 1.0}), std::invalid_argument);
    CHECK(st::makeSummaries(kLevels).intervals.levels == kLevels);
}

TEST_CASE("ensemble summary of hand-made entries", "[stats][calibration]") {
    st::PosteriorSummaries s = st::makeSummaries(kLevels);
    // errors +1 and -3, std 2 and 1, intervals of half-width 2.
    push(s, 10.0, 11.0, 2.0, 2.0);
    push(s, 20.0, 17.0, 1.0, 2.0);
    const st::EnsembleSummary e = st::summariseEnsemble(s);
    CHECK(e.entries == 2);
    CHECK_THAT(e.bias_mean, WithinULP(-1.0, 2));
    // sample std of {+1, -3} = sqrt(8), over sqrt(2) entries = 2.
    CHECK_THAT(e.bias_error, WithinRel(2.0, 1e-15));
    CHECK_THAT(e.rms_error, WithinRel(std::sqrt(5.0), 1e-15));
    CHECK_THAT(e.mean_std, WithinULP(1.5, 2));
    // sample std of {2, 1} = sqrt(0.5), over sqrt(2) entries = 0.5.
    CHECK_THAT(e.std_error, WithinRel(0.5, 1e-15));
    // pulls +0.5 and -3: mean -1.25, population std 1.75.
    CHECK_THAT(e.pull_mean, WithinULP(-1.25, 2));
    CHECK_THAT(e.pull_width, WithinRel(1.75, 1e-15));
    REQUIRE(e.coverage.size() == kLevels.size());
    CHECK_THAT(e.coverage[0], WithinULP(0.5, 0)); // only the first truth is inside

    st::PosteriorSummaries one = st::makeSummaries(kLevels);
    push(one, 10.0, 11.0, 2.0, 2.0);
    CHECK_THROWS_AS(st::summariseEnsemble(one), std::invalid_argument);
}

TEST_CASE("prior shift of hand-made entries", "[stats][calibration]") {
    st::PosteriorSummaries ref = st::makeSummaries(kLevels);
    st::PosteriorSummaries alt = st::makeSummaries(kLevels);
    push(ref, 10.0, 11.0, 2.0, 1.0);
    push(ref, 20.0, 19.0, 0.5, 1.0);
    push(alt, 10.0, 12.0, 4.0, 1.0); // shift +0.5, std ratio 2
    push(alt, 20.0, 18.0, 0.5, 1.0); // shift -2, std ratio 1
    const st::PriorShift shift = st::priorShift(ref, alt);
    CHECK_THAT(shift.shift_mean, WithinULP(-0.75, 2));
    // sample std of {+0.5, -2} = sqrt(3.125), over sqrt(2) entries = 1.25.
    CHECK_THAT(shift.shift_error, WithinRel(1.25, 1e-15));
    CHECK_THAT(shift.shift_max_abs, WithinULP(2.0, 0));
    CHECK_THAT(shift.std_ratio_mean, WithinULP(1.5, 2));

    st::PosteriorSummaries other = st::makeSummaries(kLevels);
    push(other, 10.0, 12.0, 4.0, 1.0);
    push(other, 21.0, 18.0, 0.5, 1.0);
    CHECK_THROWS_AS(st::priorShift(ref, other), std::invalid_argument);
    CHECK_THROWS_AS(st::priorShift(ref, st::makeSummaries(kLevels)), std::invalid_argument);
    st::PosteriorSummaries single_ref = st::makeSummaries(kLevels);
    st::PosteriorSummaries single_alt = st::makeSummaries(kLevels);
    push(single_ref, 10.0, 11.0, 2.0, 1.0);
    push(single_alt, 10.0, 12.0, 4.0, 1.0);
    CHECK_THROWS_AS(st::priorShift(single_ref, single_alt), std::invalid_argument);
}

TEST_CASE("pseudo-experiment event ids are disjoint blocks", "[stats][calibration]") {
    std::set<std::uint64_t> firsts;
    for (std::uint64_t section = 0; section < 3; ++section) {
        for (std::uint64_t index = 0; index < 3; ++index) {
            firsts.insert(st::experimentFirstEventId(section, index));
        }
    }
    CHECK(firsts.size() == 9);
    CHECK(st::experimentFirstEventId(0, 1) - st::experimentFirstEventId(0, 0) ==
          st::kMaxExperimentEvents + 1);
    CHECK_THROWS_AS(st::experimentFirstEventId(0, std::uint64_t{1} << st::kExperimentIndexBits),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::experimentFirstEventId(std::uint64_t{1} << 20U, 0), std::invalid_argument);
}

TEST_CASE("pseudo-experiments are reproducible and draw the truth from the domain",
          "[stats][calibration]") {
    const st::EnergyDomain narrow{11.0, 19.0, st::kEnergyGridPoints};
    const auto a = st::pseudoExperimentFromDomain(narrow, 20, kToyResolution, kSeed, 3, 7);
    const auto b = st::pseudoExperimentFromDomain(narrow, 20, kToyResolution, kSeed, 3, 7);
    const auto c = st::pseudoExperimentFromDomain(narrow, 20, kToyResolution, kSeed, 3, 8);
    CHECK_THAT(a.truth_En, WithinULP(b.truth_En, 0));
    CHECK(a.sample.reco_Ep == b.sample.reco_Ep);
    CHECK(a.truth_En >= narrow.En_min);
    CHECK(a.truth_En <= narrow.En_max);
    CHECK_FALSE(a.sample.reco_Ep == c.sample.reco_Ep);
    for (const double En : a.sample.truth_En) {
        CHECK_THAT(En, WithinULP(a.truth_En, 0));
    }
    const auto fixed = st::pseudoExperimentAtEnergy(15.0, 5, kToyResolution, kSeed, 1, 0);
    CHECK_THAT(fixed.truth_En, WithinULP(15.0, 0));
    CHECK(fixed.sample.size() == 5);
    CHECK_THROWS_AS(st::pseudoExperimentAtEnergy(15.0, st::kMaxExperimentEvents + 1, kToyResolution,
                                                 kSeed, 1, 0),
                    std::invalid_argument);
}

TEST_CASE("single-event summaries match the per-event posterior across blocks",
          "[stats][calibration]") {
    // 1100 events cross the 1024-event block boundary of the implementation.
    constexpr std::size_t kEvents = 1100;
    const Setup setup;
    const auto sample = st::pseudoSampleFromPrior(setup.domain, kEvents, kToyResolution, kSeed);
    st::PosteriorSummaries s = st::makeSummaries(kLevels);
    st::appendSingleEventSummaries(setup.model, setup.log_prior, sample, s);
    REQUIRE(s.size() == kEvents);
    REQUIRE(s.intervals.events() == kEvents);
    for (const std::size_t k :
         {std::size_t{0}, std::size_t{1023}, std::size_t{1024}, std::size_t{kEvents - 1}}) {
        const std::vector<double> log_post = st::singleEventLogPosterior(
            setup.model, sample.reco_Ep[k], sample.reco_theta_p[k], setup.log_prior);
        const st::PosteriorMoments m = st::posteriorMeanStd(setup.En, log_post);
        const auto intervals = st::credibleIntervals(setup.En, log_post.data(), kLevels);
        CHECK_THAT(s.truth[k], WithinULP(sample.truth_En[k], 0));
        CHECK_THAT(s.mean[k], WithinULP(m.mean, 0));
        CHECK_THAT(s.std[k], WithinULP(m.std, 0));
        for (std::size_t l = 0; l < kLevels.size(); ++l) {
            CHECK_THAT(s.intervals.at(k, l).lo, WithinULP(intervals[l].lo, 0));
            CHECK_THAT(s.intervals.at(k, l).hi, WithinULP(intervals[l].hi, 0));
        }
    }
}

TEST_CASE("combined summary counts the prior once", "[stats][calibration]") {
    const Setup setup;
    const auto one = st::pseudoExperimentAtEnergy(15.0, 1, kToyResolution, kSeed, 0, 0);
    st::PosteriorSummaries single = st::makeSummaries(kLevels);
    st::PosteriorSummaries combined = st::makeSummaries(kLevels);
    st::appendSingleEventSummaries(setup.model, setup.log_prior, one.sample, single);
    st::appendCombinedSummary(setup.model, setup.log_prior, one.sample, one.truth_En, combined);
    CHECK_THAT(combined.mean[0], WithinULP(single.mean[0], 0));
    CHECK_THAT(combined.std[0], WithinULP(single.std[0], 0));

    // With 100 events the combined posterior contracts around the truth.
    const auto many = st::pseudoExperimentAtEnergy(15.0, 100, kToyResolution, kSeed, 0, 1);
    st::appendCombinedSummary(setup.model, setup.log_prior, many.sample, many.truth_En, combined);
    REQUIRE(combined.size() == 2);
    CHECK(combined.std[1] < 0.2 * combined.std[0]);
    CHECK_THAT(combined.mean[1], WithinAbs(15.0, 5.0 * combined.std[1]));

    CHECK_THROWS_AS(
        st::appendCombinedSummary(setup.model, setup.log_prior, st::PseudoSample{}, 15.0, combined),
        std::invalid_argument);
    st::PosteriorSummaries no_levels;
    CHECK_THROWS_AS(
        st::appendCombinedSummary(setup.model, setup.log_prior, one.sample, 15.0, no_levels),
        std::invalid_argument);
}
