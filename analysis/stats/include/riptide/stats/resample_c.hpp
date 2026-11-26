#pragma once

// Case C experiments drawn from a Geant4 pool (D48). The pool has the flat
// source spectrum shaped by the cross section and the selection, g(E_n);
// the model of D46 describes the detected spectrum as a truncated Gaussian
// f. Each pool event is kept with probability f(E) / (c g_hat(E)), with
// g_hat the histogram of the pool truth E_n and c = max f / g_hat, so the
// kept events follow f exactly as far as g_hat follows g. The pool is
// shuffled once per row and the experiments walk it in order: an event
// visited, kept or not, is never visited again in that row, so the
// experiments of a row are disjoint and independent.

#include <cstddef>
#include <vector>

#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_c.hpp"

namespace riptide::stats {

// Bin width of g_hat. The smallest pool (+Z, ~2e5 events on 95 bins) has
// ~2000 events per bin, a ~2% error on each weight, far below what a
// Gaussian at least 0.8 MeV wide can see; finer bins would raise it.
inline constexpr double kResampleBinWidth = 0.1 * units::MeV;

// A density constant on equal bins that tile the domain, in 1/MeV.
struct BinnedDensity {
    EnergyDomain domain;
    double step = 0.0;           // MeV
    std::vector<double> density; // 1/MeV, one per bin, integrates to 1

    // Bin holding En; En_max falls in the last bin. Throws
    // std::invalid_argument outside the domain or for a non-finite En.
    [[nodiscard]] std::size_t binOf(double En) const;
};

// Histogram of `values` on bins of `step` over the domain, normalised to a
// density. Throws std::invalid_argument if the step does not tile the
// domain, a value lies outside it, or a bin is empty: an empty bin would
// make c infinite and the rejection impossible there.
[[nodiscard]] BinnedDensity histogramDensity(const std::vector<double>& values,
                                             const EnergyDomain& domain,
                                             double step = kResampleBinWidth);

// f(E) = N(E; mu, sigma) / mass on the domain, zero outside it. Throws as
// normalMassOnDomain, and for a mass below kMinTruncatedMass.
[[nodiscard]] double truncatedNormalDensity(const EnergyDomain& domain, double mu, double sigma,
                                            double En);

// c = max over the domain of f / g_hat, exact: g_hat is constant on a bin
// and f is unimodal, so on each bin the ratio peaks at the point of the bin
// closest to mu. Throws as truncatedNormalDensity.
[[nodiscard]] double resampleBound(const BinnedDensity& density, double mu, double sigma);

// The experiments of one row, as indices into the pool.
struct ResampledRow {
    std::size_t planned = 0;                           // truths given
    std::vector<std::vector<std::size_t>> experiments; // the complete ones, in order
    std::vector<std::size_t> visited;                  // pool events walked, per experiment
    std::vector<double> bounds;                        // c, per experiment
    std::size_t visited_total = 0; // including the incomplete experiment, if any

    [[nodiscard]] bool exhausted() const noexcept {
        return experiments.size() < planned;
    }
};

// One experiment of `events` kept events per truth, walking a random
// permutation of the pool drawn on `engine`; the accept draws follow on the
// same engine, so the row depends only on the pool, the truths and the
// engine. When the pool runs out the incomplete experiment is dropped and the
// later truths are not used. Throws std::invalid_argument for events == 0 or
// a pool value outside the density domain, and as resampleBound.
[[nodiscard]] ResampledRow resampleRow(const std::vector<double>& truth_En,
                                       const BinnedDensity& density,
                                       const std::vector<CaseCTruth>& truths, std::size_t events,
                                       random::Engine& engine);

} // namespace riptide::stats
