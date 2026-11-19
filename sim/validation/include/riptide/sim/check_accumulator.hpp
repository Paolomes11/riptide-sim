#pragma once

// Observed and expected quantities of the physics checks, accumulated one
// event at a time. Every expected value is computed at the event's own
// incident energy and summed, so the comparison is exact for a
// monoenergetic beam and for a flat spectrum alike, with no binning in E_n.
// No I/O, Geant4 or ROOT here: the app reads the ntuples and the evaluated
// files and hands in TruthRecords and ENDF tables, and the synthetic unit
// tests drive the same code.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "riptide/sim/beam_geometry.hpp"
#include "riptide/sim/check_stats.hpp"
#include "riptide/sim/endf_reader.hpp"
#include "riptide/sim/expected_distributions.hpp"
#include "riptide/sim/run_config.hpp"
#include "riptide/sim/truth_record.hpp"

namespace riptide::sim::checks {

// One nuclide of the medium with its evaluated cross sections.
struct Nuclide {
    Target target = Target::Hydrogen;
    std::string label;          // e.g. "H-1", "C-12", "C-0" (natural carbon)
    double number_density = 0.; // atoms per mm^3
    endf::Tab1 total;           // MT1
    endf::Tab1 elastic;         // MT2
    endf::Tab1 capture;         // MT102
};

// Reads MT1, MT2 and MT102 of `material`. Throws endf::EndfError if one is missing.
[[nodiscard]] Nuclide makeNuclide(Target target, std::string label, double number_density,
                                  const endf::EndfMaterial& material);

// The medium and the n-p angular distribution as one evaluated library sees them.
struct ReferenceModel {
    std::string library;            // e.g. "jeff33", "endfb80"
    endf::LegendreTable np_angular; // MF4/MT2 of H-1, neutron CM cosine
    std::vector<Nuclide> nuclides;

    [[nodiscard]] std::vector<double> coefficientsAt(double En) const;
    [[nodiscard]] std::vector<MediumComponent> mediumAt(double En) const;
};

// Default number of bins of the u and depth histograms. Twenty bins keep the
// expected count per bin above kMinExpectedCount already for the ~10^3-event
// smoke run (a few hundred n-p events), so even that run is compared without
// merging; production runs ask for finer bins from the command line.
inline constexpr std::size_t kDefaultBins = 20;

// Tolerance [mm] on a vertex lying outside the scintillator. It is far above
// the Geant4 surface tolerance (1e-9 mm) and the rounding of coordinates of
// order 100 mm, and far below any physical length of the checks.
inline constexpr double kVertexTolerance = 1.0e-6;

// Relative tolerance on E_n against the source range and absolute tolerance
// on each component of omega_n against the beam direction: both are copied
// through Geant4 unchanged, so only double rounding may separate them.
inline constexpr double kEnergyRelTolerance = 1.0e-12;
inline constexpr double kDirectionTolerance = 1.0e-12;

// Tolerance on u = (1 - cos_cm) / 2 above 1. With the HP target mass u <= 1
// holds exactly and rounding of the square roots stays near 1e-15; a wrong
// target mass (the bare proton instead of the H atom) overshoots by ~5e-4.
// 1e-9 lies far from both, so the counter flags only a modelling error.
inline constexpr double kUnitTolerance = 1.0e-9;

// Elastic model selection of QGSP_BIC_HP for neutrons (Geant4 11.4.1).
// G4HadronElasticPhysicsHP registers G4ParticleHPElastic up to 20 MeV
// (G4ParticleHPElastic.cc) and G4ChipsElasticModel from 19.5 MeV
// (G4HadronElasticPhysicsHP.cc). Where two models overlap,
// G4EnergyRangeManager::GetHadronicInteraction picks the one starting higher
// with probability (E - its minimum) / (overlap width), a linear ramp. The
// cross sections stay the HP ones up to 20 MeV (the last dataset added wins),
// so only the final state of n-p scattering changes in the window. CHIPS
// solves the kinematics with the bare proton mass and follows its own angular
// parametrisation, not the evaluated file.
inline constexpr double kChipsElasticMinEnergy = 19.5; // [MeV]
inline constexpr double kHpElasticMaxEnergy = 20.0;    // [MeV]

// Probability that an n-p elastic scattering at En is handed to CHIPS
// instead of the HP model: 0 below the window, 1 above it.
[[nodiscard]] double nonHpElasticProbability(double En) noexcept;

// Tolerance, relative to En, for telling the two kinematics apart. The HP
// final state reproduces hpProtonEnergy to ~1e-14 MeV and CHIPS reproduces
// the physical protonEnergy to the rounding of four-vectors of ~1 GeV, while
// the two differ by ~5.4e-4 Ep. 1e-11 En sits orders of magnitude above the
// first and below the second except for Ep below ~1e-6 MeV, where an event
// matches both and is counted as ambiguous.
inline constexpr double kModelRelTolerance = 1.0e-11;

// Number of E_n bins of the model-mix comparison across the overlap window.
inline constexpr std::size_t kModelMixBins = 10;

// Beam, box and binning of one run.
struct CheckSetup {
    Vec3 position;     // gun position [mm]
    Vec3 direction;    // unit beam direction
    Vec3 half_size;    // scintillator half edges [mm]
    double En_min = 0; // [MeV]
    double En_max = 0; // [MeV]
    BeamChord chord;   // beam path through the box; chord.length() is L
    // Rest energy [MeV] of the H-1 target in the HP elastic model (see
    // hp_recoil.hpp): AWR of the G4NDL file times the neutron mass.
    double hp_target_mass = 0;
    std::size_t bins_u = kDefaultBins;
    std::size_t bins_z = kDefaultBins;
};

// Throws std::invalid_argument if the beam misses the scintillator, a bin
// count is zero or the target mass is not positive, so a run that cannot be
// checked is never reported as passing.
[[nodiscard]] CheckSetup makeCheckSetup(const RunConfig& config, double hp_target_mass,
                                        std::size_t bins_u, std::size_t bins_z);

// Event counters. The first group can only be nonzero through a bug in the
// simulation or in the truth bookkeeping. Above the thermal threshold the HP
// target is at rest, so u > 1 and theta_p > pi/2 are kinematically
// impossible and belong there too. The second group counts known physics
// that the checks treat separately and that is reported, not cut: n-d
// elastic scattering is classified NpElastic (the truth tests Z = 1) but has
// no recoil proton, so it shows up as a first NpElastic without valid_np at
// a rate of order n_H2 / n_H1.
//
// Each n-p event with a valid Ep is classified by its kinematics: HP if Ep
// matches hpProtonEnergy, non-HP if it matches the physical protonEnergy
// (CHIPS, see kChipsElasticMinEnergy), ambiguous if it matches both (Ep near
// zero, treated as HP). An event matching neither, or a non-HP event below
// the overlap window, cannot come from the physics list as configured.
struct ConsistencyCounters {
    std::uint64_t En_outside_source = 0;
    std::uint64_t omega_not_beam = 0;
    std::uint64_t valid_np_without_scatter = 0; // valid_np with n_np_scatters < 1
    std::uint64_t vertex_outside_box = 0;       // r0 of an interacting event
    std::uint64_t invalid_proton_energy = 0;    // valid_np with Ep <= 0 or not finite
    std::uint64_t u_overflow = 0;               // u > 1 + kUnitTolerance
    std::uint64_t theta_p_beyond_half_pi = 0;
    std::uint64_t np_kinematics_unmatched = 0;  // Ep matches neither kinematics
    std::uint64_t np_non_hp_outside_window = 0; // non-HP at En <= kChipsElasticMinEnergy

    std::uint64_t np_first_without_proton = 0;
    std::uint64_t np_model_ambiguous = 0;

    [[nodiscard]] std::uint64_t mustBeZeroTotal() const {
        return En_outside_source + omega_not_beam + valid_np_without_scatter + vertex_outside_box +
               invalid_proton_energy + u_overflow + theta_p_beyond_half_pi +
               np_kinematics_unmatched + np_non_hp_outside_window;
    }
};

// Comparison of one histogram with the expectation of one library. The
// expected counts are kept per original bin; `merged` and `comparison` are
// present unless the check was skipped, and then `skipped` says why.
struct HistogramCheck {
    std::string library;
    std::vector<double> expected;
    std::optional<MergedBins> merged;
    std::optional<HistogramComparison> comparison;
    std::string skipped;
};

struct HistogramResult {
    std::vector<double> edges;
    std::vector<double> observed;
    std::vector<HistogramCheck> per_library;
};

struct FractionCheck {
    std::string library;
    std::array<double, kFirstInteractionCategories> expected{}; // expected counts
    // FirstInteraction codes entering the chi^2, in order; the categories with
    // too few expected events are folded into `absorbing_category`.
    std::vector<std::int32_t> categories;
    std::int32_t absorbing_category = 0;
    std::optional<FractionComparison> comparison;
    std::string skipped;
};

struct FractionResult {
    std::array<std::uint64_t, kFirstInteractionCategories> observed{};
    std::vector<FractionCheck> per_library;
};

struct LegendreCheck {
    std::string library;
    std::vector<double> expected; // mean of a_l(E_n) over the sample, l = 1 ... max degree
    std::vector<double> pulls;    // (<P_l> - expected) / error
};

struct LegendreResult {
    std::optional<LegendreMoments> moments;
    std::vector<LegendreCheck> per_library;
    std::string skipped;
};

struct ResidualSummary {
    std::size_t samples = 0;
    double mean = 0.0;    // [MeV]
    double rms = 0.0;     // [MeV]
    double max_abs = 0.0; // [MeV]
};

// Ep of the simulation against the physical (bare-proton) kinematics: the
// size of the HP target-mass convention, reported for the record.
struct OffsetSummary {
    ResidualSummary absolute;
    double relative_mean = 0.0; // mean of offset / physical Ep
    double relative_min = 0.0;
    double relative_max = 0.0;
};

// Share of n-p events handed to CHIPS against the ramp of
// nonHpElasticProbability. Only events classified HP or non-HP count; the
// model is drawn from E_n alone, so each is a Bernoulli trial with p(E_n):
// the total is compared through its binomial variance and, per E_n bin of
// the overlap window, through chi^2 = sum (o - e)^2 / var with ndf = bins.
struct ModelMixResult {
    std::uint64_t classified = 0;
    std::uint64_t observed_non_hp = 0;
    double expected_non_hp = 0.0;
    double variance = 0.0;     // sum of p (1 - p)
    double pull = 0.0;         // meaningful unless `skipped` is set
    std::vector<double> edges; // E_n bins [MeV] over the window inside the source range
    std::vector<double> classified_per_bin;
    std::vector<double> observed_per_bin;
    std::vector<double> expected_per_bin;
    std::vector<double> variance_per_bin;
    std::optional<double> chi2;
    int ndf = 0;
    double p_value = 0.0;
    std::string skipped;          // why the pull was not computed
    std::string chi2_skipped;     // why the per-bin chi^2 was not computed
    ResidualSummary Ep_residuals; // non-HP events against np_kinematics: rounding only
};

struct CheckResults {
    std::uint64_t events = 0;
    std::uint64_t np_events = 0;    // first interaction NpElastic with a valid recoil proton
    std::uint64_t np_hp_events = 0; // of which classified HP or ambiguous: the compared sample
    HistogramResult spectrum_u;
    HistogramResult depth_z;
    FractionResult fractions;
    LegendreResult legendre;
    ResidualSummary Ep_residuals;     // against the HP kinematics: rounding only
    OffsetSummary Ep_physical_offset; // against np_kinematics
    ModelMixResult model_mix;
    ConsistencyCounters counters;
};

class CheckAccumulator {
public:
    // Throws std::invalid_argument for no models or two with the same library name.
    CheckAccumulator(CheckSetup setup, std::vector<ReferenceModel> models);

    void add(const TruthRecord& record);

    [[nodiscard]] const CheckSetup& setup() const {
        return setup_;
    }
    [[nodiscard]] const std::vector<ReferenceModel>& models() const {
        return models_;
    }
    // Per HP n-p event, in the order added: mu = hpCosCmFromAngle(En, theta_p),
    // Ep - hpProtonEnergy(En, theta_p) and Ep - protonEnergy(En, theta_p) [MeV].
    [[nodiscard]] const std::vector<double>& cosines() const {
        return mu_;
    }
    [[nodiscard]] const std::vector<double>& EpResiduals() const {
        return Ep_residuals_;
    }
    [[nodiscard]] const std::vector<double>& EpPhysicalOffsets() const {
        return Ep_physical_offsets_;
    }

    // Runs the statistics on what was accumulated. A check without enough
    // events is reported as skipped with the reason instead of throwing.
    [[nodiscard]] CheckResults finish() const;

private:
    void addNpEvent(const TruthRecord& record);
    void addHpEvent(const TruthRecord& record);
    void addModelMix(double En, bool non_hp);
    void addDepth(const TruthRecord& record, const std::vector<double>& sigma_totals);
    [[nodiscard]] ModelMixResult finishModelMix() const;

    CheckSetup setup_;
    std::vector<ReferenceModel> models_;
    std::vector<double> u_edges_;
    std::vector<double> z_edges_;

    std::uint64_t events_ = 0;
    std::array<std::uint64_t, kFirstInteractionCategories> category_counts_{};
    std::vector<std::array<double, kFirstInteractionCategories>> category_expected_;

    std::vector<double> u_observed_;
    std::vector<std::vector<double>> u_expected_;
    std::vector<double> z_observed_;
    std::vector<std::vector<double>> z_expected_;

    std::vector<double> mu_;
    std::vector<double> Ep_residuals_;
    std::vector<double> Ep_physical_offsets_;
    std::vector<double> Ep_physical_relative_;
    std::vector<std::vector<double>> coefficient_sums_; // per library, a_l summed over n-p events
    std::uint64_t np_events_ = 0;
    std::vector<double> Ep_non_hp_residuals_;
    std::vector<double> mix_edges_; // empty when the source misses the window
    ModelMixResult mix_;            // running sums; finishModelMix adds the statistics
    ConsistencyCounters counters_;
};

} // namespace riptide::sim::checks
