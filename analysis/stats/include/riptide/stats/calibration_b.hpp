#pragma once

// Pseudo-experiment calibration of Case B (D43, D44): the per-experiment
// results of the shared-direction posterior, their ensemble numbers (error,
// declared resolution, pulls, HPD coverage, contraction with N) and the shift
// of the estimate under an alternative prior. Experiments are independent and
// each draws from its own engines, so they may run on several threads; the
// results are written by index and do not depend on the number of threads.
// Pure computation: no files, no clock.

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_b_posterior.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"

namespace riptide::stats {

// Runs task(k) for every k in [0, count), on `threads` threads: thread t takes
// k = t, t + threads, ... in increasing order (one thread runs inline, with no
// std::thread). A task must touch only its own output slot. If tasks throw,
// the exception of the lowest failing index is rethrown after every thread
// has joined, whatever the number of threads. Throws std::invalid_argument
// for threads == 0.
void parallelForEach(std::size_t count, unsigned threads,
                     const std::function<void(std::size_t)>& task);

struct CaseBExperimentResult {
    geometry::Direction omega_true;
    geometry::Direction omega_hat;
    DirectionSummary summary;
};

// The posterior of one experiment on `coarse_grid` with `log_prior`, refined
// on the default cap, and its summaries against the experiment's truth.
[[nodiscard]] CaseBExperimentResult
analyseCaseBExperiment(const CaseBModel& model, const CaseBPseudoExperiment& experiment,
                       const std::vector<geometry::Direction>& coarse_grid,
                       const DirectionLogPrior& log_prior);

// One row of the checklist: `experiments` experiments of `events` events
// each, truth isotropic, E_n uniform on truth_domain, experiment k on the ids
// of (section, k); uniform prior on the direction.
struct CaseBRowSpec {
    std::size_t events = 0;
    std::size_t experiments = 0;
    std::uint64_t section = 0;
};

[[nodiscard]] std::vector<CaseBExperimentResult>
runCaseBRow(const CaseBModel& model, const EnergyDomain& truth_domain, const CaseBRowSpec& row,
            const CaseBResolution& resolution, std::uint64_t seed,
            const std::vector<geometry::Direction>& coarse_grid, unsigned threads);

// One row on a finite pool of events with a single true direction, as a
// Geant4 sample or its toy twin (D45): experiment k is the events of
// pool at experiments[k] (from disjointExperiments), with truth truth_omega
// and the uniform direction prior. The pool is read, never copied whole.
[[nodiscard]] std::vector<CaseBExperimentResult>
runCaseBPoolRow(const CaseBModel& model, const CaseBPseudoSample& pool,
                const geometry::Direction& truth_omega,
                const std::vector<std::vector<std::size_t>>& experiments,
                const std::vector<geometry::Direction>& coarse_grid, unsigned threads);

// The alternative priors of the robustness test, on the very experiments of a
// checklist row: with `model`, a vMF prior on the direction of concentration
// kappa, its axis at angle axis_angle from the truth and a uniform azimuth
// about it; with `energy_model` (another E_n prior), the uniform direction
// prior. The azimuth of experiment k comes from the first id of
// (axis_section, k), an engine no event or truth uses.
struct CaseBPriorSpec {
    double kappa = 0.0;
    double axis_angle = 0.0; // rad
    std::uint64_t axis_section = 0;
};

struct CaseBPriorRow {
    std::vector<CaseBExperimentResult> vmf;
    std::vector<CaseBExperimentResult> energy;
};

// Throws std::invalid_argument for a negative or non-finite kappa (as
// vmfDirectionLogPrior) or an axis angle outside [0, pi].
[[nodiscard]] CaseBPriorRow
runCaseBPriorRow(const CaseBModel& model, const CaseBModel& energy_model,
                 const EnergyDomain& truth_domain, const CaseBRowSpec& row,
                 const CaseBPriorSpec& prior, const CaseBResolution& resolution, std::uint64_t seed,
                 const std::vector<geometry::Direction>& coarse_grid, unsigned threads);

// sqrt(2) * pull is d / s for a calibrated isotropic 2D Gaussian of std s per
// axis, a Rayleigh variable of unit scale: median sqrt(2 ln 2) ~ 1.177, and
// P(> 3) = exp(-4.5) ~ 0.011.
inline constexpr double kScaledPullTail = 3.0;

// The frequentist behaviour of a set of experiments.
struct CaseBEnsemble {
    std::size_t entries = 0;
    double error_mean = 0.0;
    double error_mean_error = 0.0; // sample std of the errors / sqrt(entries)
    double error_rms = 0.0;
    double sigma_mean = 0.0;       // mean declared resolution
    double sigma_mean_error = 0.0; // sample std of the sigmas / sqrt(entries)
    double pull_rms = 0.0;
    double scaled_pull_median = 0.0;        // median of sqrt(2) pull
    double scaled_pull_tail_fraction = 0.0; // fraction of sqrt(2) pull > kScaledPullTail
    double truth_in_cap_fraction = 0.0;
    std::array<double, kDirectionCredibleLevels.size()> coverage{};
};

// Needs at least two entries (for the errors of the means); else
// std::invalid_argument.
[[nodiscard]] CaseBEnsemble summariseCaseBEnsemble(const std::vector<CaseBExperimentResult>& r);

// Least-squares slope of log(rms) against log(N), as numpy.polyfit of degree
// one. Throws std::invalid_argument for sizes that differ, fewer than two
// points, a non-positive value or all N equal.
[[nodiscard]] double logLogSlope(const std::vector<double>& events, const std::vector<double>& rms);

// Coordinates of `estimate` in the tangent plane at `truth` (logarithmic map
// of the sphere), on the local basis (e_theta, e_phi) of the polar angles of
// the truth; the norm is the angle between the two. rad.
[[nodiscard]] std::array<double, 2> tangentCoordinates(const geometry::Direction& estimate,
                                                       const geometry::Direction& truth);

// The vector bias of a set of experiments (7_predictions.md, v2 note to
// section 4): the mean tangent coordinates (e_theta, e_phi) of omega_hat at
// each experiment's truth, in rad, with errors sample std / sqrt(entries).
// With a single true direction the basis is the same for every experiment,
// so the mean is a systematic offset of the estimate. Needs at least two
// entries; else std::invalid_argument.
struct DirectionBias {
    std::array<double, 2> mean{};  // rad
    std::array<double, 2> error{}; // rad
};

[[nodiscard]] DirectionBias directionBias(const std::vector<CaseBExperimentResult>& results);

// How far an alternative prior moves the estimate. The scalar shift of an
// experiment is angle(hat_ref, hat_alt) / sigma_ref; the vector shift is the
// difference of the tangent coordinates at the truth, averaged over the
// experiments and divided by the mean reference sigma. A nuisance prior that
// is not forgotten gives a large scalar shift with a vector mean near zero.
struct DirectionPriorShift {
    double shift_mean = 0.0;
    double shift_error = 0.0; // sample std of the shifts / sqrt(entries)
    double shift_max = 0.0;
    std::array<double, 2> vector_mean{};  // (e_theta, e_phi), units of mean sigma_ref
    std::array<double, 2> vector_error{}; // sample std / sqrt(entries), same units
};

// Both sets must describe the same experiments in the same order (same truth
// within kSameTruthAngle), at least two; else std::invalid_argument.
inline constexpr double kSameTruthAngle = 1e-12; // rad
[[nodiscard]] DirectionPriorShift
directionPriorShift(const std::vector<CaseBExperimentResult>& reference,
                    const std::vector<CaseBExperimentResult>& alternative);

} // namespace riptide::stats
