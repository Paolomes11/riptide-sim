#pragma once

// Output of riptide_sim_check: the results as JSON (the summary that goes into
// checks.json) and the per-bin tables as CSV. Pure formatting, no file system
// access except the atomic write, so the layout is unit-tested in every preset.
// Observed and expected values are always written side by side and unmerged;
// the merged bins appear only in the JSON, next to the chi^2 they produced.

#include <nlohmann/json.hpp>

#include <cstdint>
#include <filesystem>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/deposit_checks.hpp"
#include "riptide/sim/pstar_checks.hpp"

namespace riptide::sim::checks {

inline constexpr const char* kChecksFileName = "checks.json";
inline constexpr const char* kSpectrumUFileName = "spectrum_u.csv";
inline constexpr const char* kDepthZFileName = "depth_z.csv";
inline constexpr const char* kFractionsFileName = "fractions.csv";
inline constexpr const char* kLegendreFileName = "legendre.csv";
inline constexpr const char* kResidualsFileName = "residuals_Ep.csv";
inline constexpr const char* kModelMixFileName = "model_mix.csv";
inline constexpr const char* kDepositSpectraFileName = "deposit_spectra.csv";
inline constexpr const char* kProtonExitsFileName = "proton_exits.csv";
inline constexpr const char* kPstarComparisonFileName = "pstar_comparison.csv";

// Version of the checks.json layout and of how its numbers are formed; bumped
// on any change of its keys or of their meaning. Version 6 added
// results.deposits, version 7 results.deposits.pstar, version 8
// results.deposits.protons_nuclear_recoil (and those protons left the balance
// and PSTAR samples). Version 9 reads the rows in ascending event_id, so every
// sum, and every output byte, no longer depends on how the run's events were
// spread over the worker threads. Version 10 traces each event's own ray,
// for the far-field source of D27: depth_z is binned in the fraction of the
// event's chord, results.geometric_acceptance and the counters
// start_not_on_source and interaction_off_chord are new, and the beam object
// describes the source. On a pencil run every other number is that of 9.
inline constexpr int kChecksVersion = 10;

class ReportError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Name of a FirstInteraction code in the reports, e.g. "np_elastic".
// Throws std::out_of_range for a code outside the schema.
[[nodiscard]] std::string_view categoryName(std::int32_t code);

// The "results" object of checks.json: event counts, the five checks per
// library, the HP/CHIPS model mix and the consistency counters, the must-be-zero ones apart from
// the expected-rare ones. A skipped check carries its reason in "skipped" and no chi^2; a compared
// one has "skipped": null.
[[nodiscard]] nlohmann::json resultsJson(const CheckResults& results);

// The "beam" object of checks.json: source mode, energy range and binning;
// for a pencil beam its chord L_mm, for a far-field source the direction (or
// "isotropic"), the start disk, the expected acceptance and mean chord.
[[nodiscard]] nlohmann::json beamJson(const CheckSetup& setup);

// The "deposits" object of checks.json for a run with deposit columns: the
// setup, the event classes, the proton counts, the energy balance of the
// contained protons, the kB = 0 oracle and the counters, must-be-zero apart
// from expected-rare. "skipped" is null here; depositsSkippedJson gives the
// object of a run without deposit columns.
[[nodiscard]] nlohmann::json depositsJson(const DepositSetup& setup, const DepositResults& results);
[[nodiscard]] nlohmann::json depositsSkippedJson(const std::string& reason);

// The comparison with PSTAR: the compared and below-table protons, the
// range (mm) and light (MeV) residuals over all bins and per bin of truth_Ep.
// Undefined statistics (fewer than two protons) are null. The caller adds the
// source of the table.
[[nodiscard]] nlohmann::json pstarJson(const PstarResults& results);

// low,high,Edep,Evis: events per bin of det_Edep and det_Evis; the overflow
// above En_max is only in the JSON.
void writeDepositSpectraCsv(std::ostream& out, const DepositResults& results);
// Ep_low,Ep_high,protons,exiting: recoil protons by truth_Ep and how many left
// the scintillator.
void writeProtonExitsCsv(std::ostream& out, const DepositResults& results);

// Ep_low,Ep_high,protons,length_mean,R_csda_mean,length_residual_mean,
// length_residual_sd,Evis_mean,L_mean,Evis_residual_mean,Evis_residual_sd:
// one row per bin of truth_Ep; the means are empty without protons, the
// standard deviations with fewer than two.
void writePstarComparisonCsv(std::ostream& out, const PstarResults& results);

// low,high,observed,expected_<library>... one row per original bin.
void writeHistogramCsv(std::ostream& out, const HistogramResult& histogram);
// code,category,observed,expected_<library>...
void writeFractionsCsv(std::ostream& out, const FractionResult& fractions);
// l,value,error,expected_<library>,pull_<library>...; only the header when skipped.
void writeLegendreCsv(std::ostream& out, const LegendreResult& legendre);
// En_low,En_high,classified,observed_non_hp,expected_non_hp,variance, one row
// per E_n bin of the HP/CHIPS window; only the header when the source misses it.
void writeModelMixCsv(std::ostream& out, const ModelMixResult& mix);
// mu,Ep_residual_MeV,Ep_physical_offset_MeV, one row per HP n-p event in the
// order the events were read. Throws std::invalid_argument if the columns
// differ in length.
void writeResidualsCsv(std::ostream& out, const std::vector<double>& mu,
                       const std::vector<double>& Ep_residuals,
                       const std::vector<double>& Ep_physical_offsets);

// Writes `doc` to `path` through a temporary file and a rename, so a reader
// sees either no file or the complete one: the presence of checks.json marks a
// finished check. Throws ReportError if `path` exists or writing fails.
void writeJsonAtomically(const std::filesystem::path& path, const nlohmann::json& doc);

} // namespace riptide::sim::checks
