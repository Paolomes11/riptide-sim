#pragma once

// Pseudo-experiment calibration of Case A (D35). Samples drawn from the toy
// model are turned into per-event or combined posteriors, and these into the
// numbers of the checklist (bias, pulls, coverage), of the contraction with N
// and of the shift between two priors. Pure computation: no files, no clock.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "riptide/stats/case_a.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/numerics.hpp"
#include "riptide/stats/pseudo_experiment.hpp"
#include "riptide/stats/validate.hpp"

namespace riptide::stats {

// Mean, std and central credible intervals of a set of posteriors (one per
// event, or one per combined experiment), next to the truth each estimates.
struct PosteriorSummaries {
    std::vector<double> truth;
    std::vector<double> mean;
    std::vector<double> std;
    IntervalTable intervals; // intervals.levels fixed before the first append

    [[nodiscard]] std::size_t size() const noexcept {
        return truth.size();
    }
};

// Empty summaries for the given credible levels; throws std::invalid_argument
// for no level or a level outside (0, 1).
[[nodiscard]] PosteriorSummaries makeSummaries(const std::vector<double>& levels);

// Appends the single-event posterior (log-likelihood + log prior) of every
// event of `sample`; the truth of event k is sample.truth_En[k].
void appendSingleEventSummaries(const CaseAModel& model, const std::vector<double>& log_prior_En,
                                const PseudoSample& sample, PosteriorSummaries& out);

// Appends one entry: the combined posterior of all the events of `sample`,
// which share the energy truth_En. The prior is counted once.
void appendCombinedSummary(const CaseAModel& model, const std::vector<double>& log_prior_En,
                           const PseudoSample& sample, double truth_En, PosteriorSummaries& out);

// The frequentist behaviour of a set of summaries, with error = mean - truth.
struct EnsembleSummary {
    std::size_t entries = 0;
    double bias_mean = 0.0;  // mean error
    double bias_error = 0.0; // sample std of the errors / sqrt(entries)
    double rms_error = 0.0;
    double mean_std = 0.0;  // mean posterior std
    double std_error = 0.0; // sample std of the posterior stds / sqrt(entries)
    double pull_mean = 0.0;
    double pull_width = 0.0; // population std of the pulls
    std::vector<double> levels;
    std::vector<double> coverage;
};

// Needs at least two entries (for bias_error); else std::invalid_argument.
[[nodiscard]] EnsembleSummary summariseEnsemble(const PosteriorSummaries& summaries);

// How far an alternative prior moves the estimate, entry by entry, in units
// of the reference posterior std: (mean_alt - mean_ref) / std_ref.
struct PriorShift {
    double shift_mean = 0.0;
    double shift_error = 0.0; // sample std of the shifts / sqrt(entries)
    double shift_max_abs = 0.0;
    double std_ratio_mean = 0.0; // mean of std_alt / std_ref
};

// Both summaries must describe the same entries (same truths, same order),
// at least two (for shift_error); else std::invalid_argument.
[[nodiscard]] PriorShift priorShift(const PosteriorSummaries& reference,
                                    const PosteriorSummaries& alternative);

// Event ids of the pseudo-experiments. Experiment `index` of section
// `section` owns a block of 2^kExperimentIdBits ids: its events use the first
// ones, its truth the last one, so no two draws share an engine.
inline constexpr unsigned kExperimentIdBits = 20;
inline constexpr unsigned kExperimentIndexBits = 24;
inline constexpr std::size_t kMaxExperimentEvents = (std::size_t{1} << kExperimentIdBits) - 1;

// Throws std::invalid_argument when section or index exceed their bits.
[[nodiscard]] std::uint64_t experimentFirstEventId(std::uint64_t section, std::uint64_t index);

struct PseudoExperiment {
    double truth_En = 0.0;
    PseudoSample sample;
};

// truth_En uniform on truth_domain, then n events at that energy; n must not
// exceed kMaxExperimentEvents.
[[nodiscard]] PseudoExperiment pseudoExperimentFromDomain(const EnergyDomain& truth_domain,
                                                          std::size_t n,
                                                          const CaseAResolution& resolution,
                                                          std::uint64_t seed, std::uint64_t section,
                                                          std::uint64_t index);

[[nodiscard]] PseudoExperiment pseudoExperimentAtEnergy(double En, std::size_t n,
                                                        const CaseAResolution& resolution,
                                                        std::uint64_t seed, std::uint64_t section,
                                                        std::uint64_t index);

} // namespace riptide::stats
