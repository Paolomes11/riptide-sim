#pragma once

// Expected outcome of the physics checks of a run, computed from the
// evaluated data alone, before the run exists. The checks accumulate their
// expectations event by event at the sampled E_n; here the same per-energy
// quantities are averaged over the source spectrum instead, so a prediction
// and a check of the same run agree up to the statistics of the sample. The
// spectrum average is a quadrature in E_n: one node for a monoenergetic
// source, the midpoint rule for a flat one, with its error measured by
// doubling the nodes (predictionShift). A far-field source gives each event
// its own chord: the per-chord probabilities are averaged as well, over the
// chord distribution of the source (ChordDistribution).
//
// Pure functions, no I/O, Geant4 or ROOT: the app hands in the reference
// models and writes the JSON.

#include <nlohmann/json.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/expected_distributions.hpp"

namespace riptide::sim::checks {

// Nodes and weights of a quadrature rule.
struct QuadratureRule {
    std::vector<double> nodes;
    std::vector<double> weights;
};

// Gauss-Legendre rule with `points` nodes on [-1, 1]: exact for polynomials
// of degree up to 2 points - 1. Throws std::invalid_argument for 0 points.
[[nodiscard]] QuadratureRule gaussLegendre(std::size_t points);

// Spectrum average over [En_min, En_max] with weights summing to 1: a single
// node at En_min when the range is empty (monoenergetic source), else the
// midpoint rule with `points` nodes (flat source). Throws
// std::invalid_argument for En_max < En_min, a non-positive En_min or 0 points.
[[nodiscard]] QuadratureRule sourceQuadrature(double En_min, double En_max, std::size_t points);

// <P_degree(mu)^2> under f(mu) = 1/2 + sum_l (2l+1)/2 a_l P_l(mu), exact by
// Gauss-Legendre since the integrand is a polynomial. With <P_degree> = a_degree
// it gives the per-event variance of P_degree, hence the error of a Legendre
// moment before the run. Throws std::invalid_argument for degree < 1.
[[nodiscard]] double legendreSecondMoment(const std::vector<double>& a, int degree);

// z such that P(|Z| > z) = p for a standard normal Z. Throws
// std::invalid_argument for p outside (0, 1).
[[nodiscard]] double normalTwoSidedQuantile(double p);

// Spectrum averages of one reference library.
struct LibraryPrediction {
    std::string library;
    double sigma_total = 0.0; // macroscopic total cross section [1/mm]
    // P(first interaction category), indexed by schema::FirstInteraction;
    // entry Unknown is P(no interaction).
    std::array<double, kFirstInteractionCategories> probabilities{};
    double hp_np_probability = 0.0;     // P(first n-p, final state by HP)
    double non_hp_np_probability = 0.0; // P(first n-p, final state by CHIPS)
    double non_hp_variance = 0.0;       // spectrum average of P(n-p) p (1 - p), p the CHIPS share
    // Over the HP n-p events, l = 1 ... max degree: mean of a_l(E_n), the
    // expected <P_l>, and the spread of P_l(mu) per event.
    std::vector<double> legendre;
    std::vector<double> legendre_spread;
    std::vector<double> u_probabilities; // per bin of u over [0, 1], HP n-p events
};

// Distribution of the chord of an event through the box. Pencil: the one
// beam chord with weight 1. Far-field: the chords of `samples` rays drawn
// from the source with `seed`, grouped in bins of kChordBinWidth, each bin
// at the mean length of its rays. The weights of the hit chords sum to the
// closed-form acceptance and miss_weight is its complement, so the sample
// only shapes the chord distribution and the miss fraction carries no
// sampling error.
struct ChordDistribution {
    bool far_field = false;
    std::vector<double> lengths; // [mm], ascending
    std::vector<double> weights;
    double miss_weight = 0.0;
    std::uint64_t samples = 0;
    std::uint64_t seed = 0;
    std::uint64_t sample_hits = 0;
    double sample_mean = 0.0; // mean chord of the sampled hits [mm]
    double sample_sd = 0.0;   // their sample standard deviation [mm]
    double acceptance = 1.0;  // closed form, as in CheckSetup
    double expected_mean_chord = 0.0;
};

// Width [mm] of the chord bins. Inside a bin exp(-Sigma L) is replaced by its
// value at the bin mean, a relative error of order (Sigma w)^2 / 24: far below
// 1e-8 for Sigma of a few 1e-2 / mm.
inline constexpr double kChordBinWidth = 1.0e-2;

// The pencil chord of setup. Throws std::invalid_argument for a far-field setup.
[[nodiscard]] ChordDistribution pencilChord(const CheckSetup& setup);

// Far-field chords sampled from the source of setup with a seeded
// std::mt19937_64. Throws std::invalid_argument for a pencil setup, 0
// samples or a sample with no ray crossing the box.
[[nodiscard]] ChordDistribution farFieldChords(const CheckSetup& setup, std::uint64_t samples,
                                               std::uint64_t seed);

struct RunPrediction {
    std::uint64_t events = 0;
    ChordDistribution chords;
    double En_min = 0.0;
    double En_max = 0.0;
    std::size_t energy_nodes = 0;
    bool overlaps_model_window = false; // the source reaches the HP/CHIPS overlap
    // Relative offset (hpProtonEnergy - protonEnergy) / protonEnergy over
    // theta_p in [0, pi/2) and the source range: the range the HP events of
    // the run must fall in.
    double offset_relative_min = 0.0;
    double offset_relative_max = 0.0;
    std::vector<LibraryPrediction> libraries; // in the order of the models
};

// Throws std::invalid_argument for no models, 0 events, 0 energy nodes, an
// empty chord distribution or a source without HP n-p events.
[[nodiscard]] RunPrediction predictRun(const CheckSetup& setup,
                                       const std::vector<ReferenceModel>& models,
                                       const ChordDistribution& chords, std::uint64_t events,
                                       std::size_t energy_nodes);

// The same along the pencil chord of setup (throws for a far-field setup).
[[nodiscard]] RunPrediction predictRun(const CheckSetup& setup,
                                       const std::vector<ReferenceModel>& models,
                                       std::uint64_t events, std::size_t energy_nodes);

// Statistical errors expected for the run under library `reference`.
struct ExpectedErrors {
    double hp_np_events = 0.0;
    std::array<double, kFirstInteractionCategories> fraction_sigma{}; // sqrt(p (1 - p) / N)
    std::vector<double> legendre_sigma;                               // spread / sqrt(N_hp)
};

[[nodiscard]] ExpectedErrors expectedErrors(const RunPrediction& run, std::size_t reference);

// Largest difference between two predictions of the same run, in units of
// the expected statistical errors of `reference`: used for the quadrature
// error (same models, doubled nodes) and for the separation between two
// libraries. `u_noncentrality` is the chi^2 excess sum N_hp (q - q_ref)^2 /
// q_ref that the u histogram would show if the sample followed `other`.
struct PredictionShift {
    double fractions = 0.0; // max over categories of |dp| / sigma
    double legendre = 0.0;  // max over l of |da_l| / sigma
    double u_noncentrality = 0.0;
};

// Throws std::invalid_argument if the predictions have different binning,
// event counts or degrees.
[[nodiscard]] PredictionShift predictionShift(const LibraryPrediction& reference,
                                              const LibraryPrediction& other,
                                              const ExpectedErrors& errors, std::uint64_t events);

// Statistical acceptance tests of a run against one library: the chi^2 of
// u, depth and fractions, one pull per Legendre degree, when the source
// reaches the overlap window the model-mix pull and chi^2 and, for a
// far-field source, the acceptance pull and the mean-chord pull (the latter
// only when the chords have a spread). The family-wise
// false-alarm probability is split evenly (Bonferroni), which stays
// conservative for the correlated Legendre pulls.
struct TestPlan {
    int chi2_tests = 0;
    int pull_tests = 0;
    double alpha_family = 0.0;
    double alpha_per_test = 0.0; // a chi^2 passes if its p-value exceeds it
    double pull_threshold = 0.0; // a pull passes if |pull| stays below it
};

// Throws std::invalid_argument for alpha_family outside (0, 1).
[[nodiscard]] TestPlan testPlan(const RunPrediction& run, double alpha_family);

// The whole prediction with derived counts, errors, the shift of every other
// library from the first and the test plan.
[[nodiscard]] nlohmann::json
predictionJson(const RunPrediction& run, const PredictionShift& quadrature, const TestPlan& plan);

// The "chords" object of the prediction: mode, closed-form acceptance and
// mean chord and, for a far-field source, the sample that shaped the
// distribution.
[[nodiscard]] nlohmann::json chordsJson(const ChordDistribution& chords);

} // namespace riptide::sim::checks
