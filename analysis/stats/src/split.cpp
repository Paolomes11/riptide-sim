#include "riptide/stats/split.hpp"

#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace riptide::stats {

namespace {

void requirePositiveSize(std::size_t events_per_experiment) {
    if (events_per_experiment == 0) {
        throw std::invalid_argument("split: events per experiment must be positive");
    }
}

} // namespace

std::size_t maxDisjointExperiments(std::size_t events, std::size_t events_per_experiment) {
    requirePositiveSize(events_per_experiment);
    return events / events_per_experiment;
}

std::vector<std::vector<std::size_t>> disjointExperiments(std::size_t events,
                                                          std::size_t events_per_experiment,
                                                          std::size_t experiments,
                                                          random::Engine& engine) {
    // Comparing with the quotient also rules out an overflow of the product.
    if (experiments > maxDisjointExperiments(events, events_per_experiment)) {
        throw std::invalid_argument("split: the experiments do not fit in the sample");
    }
    const std::size_t used = experiments * events_per_experiment;

    std::vector<std::size_t> order(events);
    std::iota(order.begin(), order.end(), std::size_t{0});
    // Position i takes a uniform pick among the events not yet placed; the
    // positions after `used` are never read, so they are not shuffled.
    for (std::size_t i = 0; i < used; ++i) {
        const auto offset =
            static_cast<std::size_t>(random::uniformBelow(engine, std::uint64_t{events - i}));
        std::swap(order[i], order[i + offset]);
    }

    std::vector<std::vector<std::size_t>> blocks(experiments);
    for (std::size_t m = 0; m < experiments; ++m) {
        const auto first = order.begin() + static_cast<std::ptrdiff_t>(m * events_per_experiment);
        blocks[m].assign(first, first + static_cast<std::ptrdiff_t>(events_per_experiment));
    }
    return blocks;
}

PseudoSample pseudoSampleSubset(const PseudoSample& sample,
                                const std::vector<std::size_t>& indices) {
    const std::size_t n = sample.size();
    if (sample.truth_Ep.size() != n || sample.truth_theta_p.size() != n ||
        sample.reco_Ep.size() != n || sample.reco_theta_p.size() != n) {
        throw std::invalid_argument("pseudoSampleSubset: columns of different lengths");
    }
    PseudoSample subset;
    subset.truth_En.reserve(indices.size());
    subset.truth_Ep.reserve(indices.size());
    subset.truth_theta_p.reserve(indices.size());
    subset.reco_Ep.reserve(indices.size());
    subset.reco_theta_p.reserve(indices.size());
    for (const std::size_t k : indices) {
        if (k >= n) {
            throw std::invalid_argument("pseudoSampleSubset: index out of range");
        }
        subset.truth_En.push_back(sample.truth_En[k]);
        subset.truth_Ep.push_back(sample.truth_Ep[k]);
        subset.truth_theta_p.push_back(sample.truth_theta_p[k]);
        subset.reco_Ep.push_back(sample.reco_Ep[k]);
        subset.reco_theta_p.push_back(sample.reco_theta_p[k]);
    }
    return subset;
}

CaseBPseudoSample caseBPseudoSampleSubset(const CaseBPseudoSample& sample,
                                          const std::vector<std::size_t>& indices) {
    const std::size_t n = sample.size();
    if (sample.truth_Ep.size() != n || sample.truth_theta_p.size() != n ||
        sample.truth_track.size() != n || sample.reco_Ep.size() != n ||
        sample.reco_track.size() != n) {
        throw std::invalid_argument("caseBPseudoSampleSubset: columns of different lengths");
    }
    CaseBPseudoSample subset;
    subset.truth_En.reserve(indices.size());
    subset.truth_Ep.reserve(indices.size());
    subset.truth_theta_p.reserve(indices.size());
    subset.truth_track.reserve(indices.size());
    subset.reco_Ep.reserve(indices.size());
    subset.reco_track.reserve(indices.size());
    for (const std::size_t k : indices) {
        if (k >= n) {
            throw std::invalid_argument("caseBPseudoSampleSubset: index out of range");
        }
        subset.truth_En.push_back(sample.truth_En[k]);
        subset.truth_Ep.push_back(sample.truth_Ep[k]);
        subset.truth_theta_p.push_back(sample.truth_theta_p[k]);
        subset.truth_track.push_back(sample.truth_track[k]);
        subset.reco_Ep.push_back(sample.reco_Ep[k]);
        subset.reco_track.push_back(sample.reco_track[k]);
    }
    return subset;
}

} // namespace riptide::stats
