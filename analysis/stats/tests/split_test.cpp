#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

#include "riptide/common/random.hpp"
#include "riptide/stats/split.hpp"

namespace st = riptide::stats;
using riptide::random::Engine;

namespace {

constexpr std::uint64_t kSeed = 20260907;
// Frequency checks at 4 standard errors: a correct shuffle fails one with
// probability ~6e-5, and the engine sequence is fixed, so a pass is stable.
constexpr double kMcSigmas = 4.0;

std::vector<std::size_t> flatten(const std::vector<std::vector<std::size_t>>& blocks) {
    std::vector<std::size_t> all;
    for (const auto& block : blocks) {
        all.insert(all.end(), block.begin(), block.end());
    }
    return all;
}

} // namespace

TEST_CASE("disjointExperiments gives M blocks of N distinct in-range indices", "[stats][split]") {
    constexpr std::size_t kEvents = 1000;
    constexpr std::size_t kPerExperiment = 30;
    constexpr std::size_t kExperiments = 33;
    Engine engine{kSeed};
    const auto blocks = st::disjointExperiments(kEvents, kPerExperiment, kExperiments, engine);
    REQUIRE(blocks.size() == kExperiments);
    for (const auto& block : blocks) {
        CHECK(block.size() == kPerExperiment);
    }
    const auto all = flatten(blocks);
    CHECK(std::all_of(all.begin(), all.end(), [](std::size_t k) { return k < kEvents; }));
    const std::set<std::size_t> unique(all.begin(), all.end());
    CHECK(unique.size() == all.size());
}

TEST_CASE("disjointExperiments uses every event when the blocks fill the sample",
          "[stats][split]") {
    constexpr std::size_t kPerExperiment = 7;
    constexpr std::size_t kExperiments = 11;
    constexpr std::size_t kEvents = kPerExperiment * kExperiments;
    Engine engine{kSeed};
    auto all = flatten(st::disjointExperiments(kEvents, kPerExperiment, kExperiments, engine));
    std::sort(all.begin(), all.end());
    std::vector<std::size_t> expected(kEvents);
    std::iota(expected.begin(), expected.end(), std::size_t{0});
    CHECK(all == expected);
}

TEST_CASE("disjointExperiments is reproducible and seed dependent", "[stats][split]") {
    constexpr std::size_t kEvents = 500;
    constexpr std::size_t kPerExperiment = 10;
    constexpr std::size_t kExperiments = 20;
    Engine a{kSeed};
    Engine b{kSeed};
    Engine c{kSeed + 1};
    const auto first = st::disjointExperiments(kEvents, kPerExperiment, kExperiments, a);
    CHECK(first == st::disjointExperiments(kEvents, kPerExperiment, kExperiments, b));
    CHECK(first != st::disjointExperiments(kEvents, kPerExperiment, kExperiments, c));
}

TEST_CASE("disjointExperiments rejects blocks that do not fit", "[stats][split]") {
    Engine engine{kSeed};
    CHECK_THROWS_AS(st::disjointExperiments(10, 0, 1, engine), std::invalid_argument);
    CHECK_THROWS_AS(st::disjointExperiments(10, 3, 4, engine), std::invalid_argument);
    // M * N wraps around in std::size_t: the check must not multiply.
    constexpr std::size_t kHuge = std::numeric_limits<std::size_t>::max() / 2 + 1;
    CHECK_THROWS_AS(st::disjointExperiments(10, kHuge, 2, engine), std::invalid_argument);
    CHECK(st::disjointExperiments(10, 3, 0, engine).empty());
    CHECK(st::disjointExperiments(0, 3, 0, engine).empty());
}

TEST_CASE("maxDisjointExperiments is the integer quotient", "[stats][split]") {
    CHECK(st::maxDisjointExperiments(10, 3) == 3);
    CHECK(st::maxDisjointExperiments(9, 3) == 3);
    CHECK(st::maxDisjointExperiments(2, 3) == 0);
    CHECK_THROWS_AS(st::maxDisjointExperiments(10, 0), std::invalid_argument);
}

TEST_CASE("disjointExperiments picks uniform ordered subsets", "[stats][split]") {
    // One block of N out of n: each event is in the block with probability
    // N / n and sits at a given position with probability 1 / n.
    constexpr std::size_t kEvents = 10;
    constexpr std::size_t kPerExperiment = 3;
    constexpr std::size_t kRepeats = 200000;
    Engine engine{kSeed};
    std::array<std::size_t, kEvents> in_block{};
    std::array<std::size_t, kPerExperiment> event0_at{};
    for (std::size_t r = 0; r < kRepeats; ++r) {
        const auto block = st::disjointExperiments(kEvents, kPerExperiment, 1, engine).front();
        for (std::size_t pos = 0; pos < kPerExperiment; ++pos) {
            ++in_block.at(block[pos]);
            if (block[pos] == 0) {
                ++event0_at.at(pos);
            }
        }
    }
    const double n = static_cast<double>(kRepeats);
    const double p_in = static_cast<double>(kPerExperiment) / static_cast<double>(kEvents);
    const double p_at = 1.0 / static_cast<double>(kEvents);
    for (const std::size_t count : in_block) {
        CHECK_THAT(static_cast<double>(count) / n,
                   Catch::Matchers::WithinAbs(p_in, kMcSigmas * std::sqrt(p_in * (1 - p_in) / n)));
    }
    for (const std::size_t count : event0_at) {
        CHECK_THAT(static_cast<double>(count) / n,
                   Catch::Matchers::WithinAbs(p_at, kMcSigmas * std::sqrt(p_at * (1 - p_at) / n)));
    }
}

TEST_CASE("pseudoSampleSubset copies the selected events in order", "[stats][split]") {
    st::PseudoSample sample;
    constexpr std::size_t kEvents = 5;
    for (std::size_t k = 0; k < kEvents; ++k) {
        const auto x = static_cast<double>(k);
        sample.truth_En.push_back(10.0 + x);
        sample.truth_Ep.push_back(1.0 + x);
        sample.truth_theta_p.push_back(0.1 * x);
        sample.reco_Ep.push_back(2.0 + x);
        sample.reco_theta_p.push_back(0.2 * x);
    }
    const std::vector<std::size_t> indices{3, 0, 4};
    const auto subset = st::pseudoSampleSubset(sample, indices);
    REQUIRE(subset.size() == indices.size());
    REQUIRE(subset.reco_theta_p.size() == indices.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
        const std::size_t k = indices[i];
        CHECK_THAT(subset.truth_En[i], Catch::Matchers::WithinULP(sample.truth_En[k], 0));
        CHECK_THAT(subset.truth_Ep[i], Catch::Matchers::WithinULP(sample.truth_Ep[k], 0));
        CHECK_THAT(subset.truth_theta_p[i], Catch::Matchers::WithinULP(sample.truth_theta_p[k], 0));
        CHECK_THAT(subset.reco_Ep[i], Catch::Matchers::WithinULP(sample.reco_Ep[k], 0));
        CHECK_THAT(subset.reco_theta_p[i], Catch::Matchers::WithinULP(sample.reco_theta_p[k], 0));
    }

    CHECK_THROWS_AS(st::pseudoSampleSubset(sample, {kEvents}), std::invalid_argument);
    sample.reco_Ep.pop_back();
    CHECK_THROWS_AS(st::pseudoSampleSubset(sample, {0}), std::invalid_argument);
}

TEST_CASE("caseBPseudoSampleSubset copies the selected events in order", "[stats][split]") {
    st::CaseBPseudoSample sample;
    constexpr std::size_t kEvents = 5;
    for (std::size_t k = 0; k < kEvents; ++k) {
        const auto x = static_cast<double>(k);
        sample.truth_En.push_back(10.0 + x);
        sample.truth_Ep.push_back(1.0 + x);
        sample.truth_theta_p.push_back(0.1 * x);
        sample.truth_track.push_back({x, 0.0, 1.0});
        sample.reco_Ep.push_back(2.0 + x);
        sample.reco_track.push_back({0.0, x, 1.0});
    }
    const std::vector<std::size_t> indices{3, 0, 4};
    const auto subset = st::caseBPseudoSampleSubset(sample, indices);
    REQUIRE(subset.size() == indices.size());
    REQUIRE(subset.reco_track.size() == indices.size());
    using Catch::Matchers::WithinULP;
    for (std::size_t i = 0; i < indices.size(); ++i) {
        const std::size_t k = indices[i];
        CHECK_THAT(subset.truth_En[i], WithinULP(sample.truth_En[k], 0));
        CHECK_THAT(subset.truth_Ep[i], WithinULP(sample.truth_Ep[k], 0));
        CHECK_THAT(subset.truth_theta_p[i], WithinULP(sample.truth_theta_p[k], 0));
        CHECK_THAT(subset.truth_track[i].x, WithinULP(sample.truth_track[k].x, 0));
        CHECK_THAT(subset.reco_Ep[i], WithinULP(sample.reco_Ep[k], 0));
        CHECK_THAT(subset.reco_track[i].y, WithinULP(sample.reco_track[k].y, 0));
    }

    CHECK_THROWS_AS(st::caseBPseudoSampleSubset(sample, {kEvents}), std::invalid_argument);
    sample.reco_track.pop_back();
    CHECK_THROWS_AS(st::caseBPseudoSampleSubset(sample, {0}), std::invalid_argument);
}
