#pragma once

// Checks of the deposit, visible-energy and recoil-proton columns of a schema
// v2 run, accumulated one event at a time. Like check_accumulator.hpp it has no
// I/O, Geant4 or ROOT: the app reads the rows and hands in the records, and the
// synthetic unit tests drive the same code. Comparisons with external tables
// of range and light yield are not here; these are the deterministic oracles
// and the distributions reported without a threshold.

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "riptide/common/event_schema.hpp"
#include "riptide/sim/deposit_record.hpp"
#include "riptide/sim/run_config.hpp"
#include "riptide/sim/truth_record.hpp"

namespace riptide::sim::checks {

inline constexpr std::size_t kEventClasses = 7; // EventClass codes 0 ... 6

// First event schema with the deposit columns; older runs skip these checks.
inline constexpr std::int32_t kFirstDepositSchemaVersion = 2;

// Tolerance [MeV] on sums of per-step deposits compared with each other: each
// side sums at most a few thousand doubles of order 1 MeV, so rounding stays
// far below it, and any physical deposit far above.
inline constexpr double kDepositSumTolerance = 1.0e-9;

// Margin [MeV] by which a recoil proton may deposit more than its creation
// energy. G4NuclearStopping (Geant4 11.4.1) caps the loss of a step at the
// pre-step kinetic energy rather than at what ionisation left, and
// G4ParticleChangeForLoss then clamps the energy to 0 while depositing the
// whole loss, so the last step of a stopping proton can overshoot by up to one
// step of nuclear stopping (1.9e-4 MeV seen at 14 MeV). 0.01 MeV is far above
// that and far below the MeV-scale error of following the wrong track; the
// excesses inside it are counted as known physics.
inline constexpr double kNuclearStoppingSlack = 0.01;

// Distance [mm] within which an end point counts as lying on the scintillator
// surface: far above the Geant4 surface tolerance (1e-9 mm) and the rounding
// of coordinates of order 100 mm, far below the micrometre-scale last step of
// a stopping proton.
inline constexpr double kSurfaceTolerance = 1.0e-6;

// Bins of the det_Edep/det_Evis spectra and of the proton exit fraction.
inline constexpr std::size_t kDefaultDepositBins = 50;

struct DepositSetup {
    Vec3 half_size; // scintillator half edges [mm], box centred at the origin
    double En_max =
        0.0; // upper edge of the source [MeV]: the spectra and the Ep bins span [0, En_max]
    double birks_kB = 0.0; // [mm/MeV]; 0 turns on the Evis == Edep oracle
    std::size_t bins = kDefaultDepositBins;
};

// Throws std::invalid_argument for a non-positive box edge or En_max, a
// negative kB or zero bins, or a config without the scintillator block
// (a version 1 configuration has no deposit columns to check).
[[nodiscard]] DepositSetup makeDepositSetup(const RunConfig& config,
                                            std::size_t bins = kDefaultDepositBins);

// The first group can only be nonzero through a bug in the simulation or in
// the bookkeeping; the second counts known physics, reported and not cut.
struct DepositCounters {
    std::uint64_t negative_deposit = 0;            // any of det_/p_ Edep, Evis below -tolerance
    std::uint64_t Evis_above_Edep = 0;             // det_Evis > det_Edep
    std::uint64_t p_Evis_above_p_Edep = 0;         // truth_p_Evis > truth_p_Edep
    std::uint64_t p_Edep_above_det_Edep = 0;       // the proton is part of the event total
    std::uint64_t p_Edep_above_Ep_slack = 0;       // beyond kNuclearStoppingSlack
    std::uint64_t proton_without_track = 0;        // valid_np with a non-positive path length
    std::uint64_t deposit_without_interaction = 0; // primary never interacted, Edep > 0
    std::uint64_t class_mismatch = 0;              // stored class differs from the counters
    std::uint64_t first_interaction_uncounted = 0; // first interaction kind has a zero counter
    std::uint64_t exit_not_on_surface = 0;         // p_exits with an end point off the surface
    std::uint64_t end_outside_box = 0;             // proton end point outside the box
    std::uint64_t kB_zero_visible_mismatch = 0;    // kB = 0 and |Evis - Edep| > tolerance

    std::uint64_t p_Edep_above_Ep = 0; // within the slack: G4NuclearStopping at end of range
    std::uint64_t contained_end_on_surface = 0; // stopped exactly on the surface

    [[nodiscard]] std::uint64_t mustBeZeroTotal() const {
        return negative_deposit + Evis_above_Edep + p_Evis_above_p_Edep + p_Edep_above_det_Edep +
               p_Edep_above_Ep_slack + proton_without_track + deposit_without_interaction +
               class_mismatch + first_interaction_uncounted + exit_not_on_surface +
               end_outside_box + kB_zero_visible_mismatch;
    }
};

// truth_p_Edep - truth_Ep of recoil protons that stop inside the scintillator
// without a hadronic interaction: the whole creation energy should be
// deposited by the track, up to the G4NuclearStopping excess (positive) and
// energy carried off by secondaries above the production cut (negative).
struct BalanceSummary {
    std::uint64_t samples = 0;
    double mean = 0.0;          // [MeV]
    double min = 0.0;           // [MeV]
    double max = 0.0;           // [MeV]
    std::uint64_t positive = 0; // above +kDepositSumTolerance
    std::uint64_t negative = 0; // below -kDepositSumTolerance
};

struct DepositResults {
    std::uint64_t events = 0;
    std::uint64_t np_events = 0; // valid_np: a recoil proton was followed
    std::array<std::uint64_t, kEventClasses> classes{};
    std::uint64_t events_secondary_n_hadronic = 0; // truth_n_secondary_n_hadronic > 0
    std::uint64_t protons_hadronic = 0;            // truth_p_n_hadronic > 0
    std::uint64_t protons_exiting = 0;
    BalanceSummary contained_balance;
    double kB_zero_max_abs_difference = 0.0; // max |det_Evis - det_Edep| [MeV], all events
    // [0, En_max] in `bins` equal bins; deposits above En_max (captures with
    // their Q-value) go to the overflow counters, not to the last bin.
    std::vector<double> edges;
    std::vector<std::uint64_t> Edep_counts; // events with det_Edep > 0
    std::vector<std::uint64_t> Evis_counts; // events with det_Evis > 0
    std::uint64_t Edep_overflow = 0;
    std::uint64_t Evis_overflow = 0;
    std::vector<std::uint64_t> protons_per_Ep; // by truth_Ep, same edges
    std::vector<std::uint64_t> exiting_per_Ep;
    DepositCounters counters;
};

class DepositAccumulator {
public:
    explicit DepositAccumulator(DepositSetup setup);

    // `stored_class` is the truth_event_class column as written, compared with
    // the class the counters of `truth` give.
    void add(const TruthRecord& truth, const DepositRecord& deposit, std::int32_t stored_class);

    [[nodiscard]] const DepositSetup& setup() const {
        return setup_;
    }
    [[nodiscard]] DepositResults finish() const;

private:
    void addProton(const TruthRecord& truth, const DepositRecord& deposit);
    void addSpectra(const DepositRecord& deposit);
    [[nodiscard]] bool onSurface(const Vec3& p) const;
    [[nodiscard]] bool insideBox(const Vec3& p) const;

    DepositSetup setup_;
    DepositResults results_;
    double balance_sum_ = 0.0;
};

} // namespace riptide::sim::checks
