#pragma once

// Comparison of the recoil protons of a schema v2 run with the NIST PSTAR
// table: the path length against the CSDA range and the visible energy
// against Birks' law integrated over the tabulated stopping power, at the
// kB of the run. Only protons that stop inside the scintillator without a
// hadronic interaction enter, the sample of the contained energy balance:
// for them the whole track is in the columns. The expectation is computed at
// each proton's own energy, so the width of the Ep bins does not enter.
//
// No thresholds here: these are the numbers the predictions judge.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "riptide/sim/deposit_checks.hpp"
#include "riptide/sim/deposit_record.hpp"
#include "riptide/sim/pstar.hpp"
#include "riptide/sim/truth_record.hpp"

namespace riptide::sim::checks {

// Observed and expected values of one quantity and their difference d =
// observed - expected, summed with Welford's update so that the variance
// stays accurate over 1e7 samples of nearly equal values.
class ResidualStats {
public:
    void add(double observed, double expected);

    [[nodiscard]] std::uint64_t samples() const {
        return samples_;
    }
    // Means; 0 without samples.
    [[nodiscard]] double observedMean() const {
        return observed_mean_;
    }
    [[nodiscard]] double expectedMean() const {
        return expected_mean_;
    }
    [[nodiscard]] double residualMean() const {
        return residual_mean_;
    }
    // Sample standard deviation of d and the standard error of its mean;
    // empty with fewer than two samples.
    [[nodiscard]] std::optional<double> residualSd() const;
    [[nodiscard]] std::optional<double> residualSem() const;
    // Mean of d over the mean expectation, sum(d) / sum(expected), and its
    // standard error; empty without samples or with a non-positive expectation
    // (the error also with fewer than two samples).
    [[nodiscard]] std::optional<double> relativeResidual() const;
    [[nodiscard]] std::optional<double> relativeSem() const;

private:
    std::uint64_t samples_ = 0;
    double observed_mean_ = 0.0;
    double expected_mean_ = 0.0;
    double residual_mean_ = 0.0;
    double residual_m2_ = 0.0; // sum of squared deviations from the running mean
};

struct PstarResults {
    std::uint64_t samples = 0;     // compared protons
    std::uint64_t below_table = 0; // contained, Ep below the first tabulated energy
    // [0, En_max] in `bins` equal bins of truth_Ep, the edges of the deposit checks.
    std::vector<double> edges;
    std::vector<ResidualStats> range_bins; // truth_p_length vs R_CSDA(Ep) [mm]
    std::vector<ResidualStats> light_bins; // truth_p_Evis vs L(Ep) [MeV]
    ResidualStats range;                   // all bins together
    ResidualStats light;
};

class PstarComparison {
public:
    // kB from the setup, in mm/MeV like the table's MeV/mm. Throws
    // std::invalid_argument for a non-positive En_max or zero bins, and
    // PstarError for an invalid table or kB, or En_max above the table.
    PstarComparison(const DepositSetup& setup, PstarTable table);

    // Selects the proton itself: valid_np, not exiting, no hadronic interaction.
    void add(const TruthRecord& truth, const DepositRecord& deposit);

    [[nodiscard]] PstarResults finish() const {
        return results_;
    }

private:
    double En_max_;
    std::size_t bins_;
    BirksLightTable light_;
    PstarResults results_;
};

} // namespace riptide::sim::checks
