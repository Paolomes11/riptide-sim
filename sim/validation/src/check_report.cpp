#include "riptide/sim/check_report.hpp"

#include <array>
#include <cstddef>
#include <fstream>
#include <ios>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>

#include "riptide/common/event_schema.hpp"

namespace riptide::sim::checks {

namespace {

constexpr int kJsonIndent = 2;
constexpr const char* kTemporarySuffix = ".tmp";

// Same order as the schema::FirstInteraction codes 0 ... 4.
constexpr std::array<std::string_view, kFirstInteractionCategories> kCategoryNames{
    "unknown", "np_elastic", "nc_elastic", "nc_inelastic", "other"};

// Same order as the schema::EventClass codes 0 ... 6.
constexpr std::array<std::string_view, kEventClasses> kEventClassNames{
    "no_interaction",     "single_np", "multiple_np", "nc_elastic_only",
    "nc_inelastic_no_np", "mixed",     "other"};

// Enough digits for every double to read back to the same value.
void setFullPrecision(std::ostream& out) {
    out.precision(std::numeric_limits<double>::max_digits10);
}

nlohmann::json skippedJson(const std::string& reason) {
    return reason.empty() ? nlohmann::json(nullptr) : nlohmann::json(reason);
}

double sum(const std::vector<double>& values) {
    double total = 0.0;
    for (const double v : values) {
        total += v;
    }
    return total;
}

nlohmann::json histogramJson(const HistogramResult& histogram) {
    nlohmann::json doc;
    doc["bins"] = histogram.observed.size();
    doc["low"] = histogram.edges.front();
    doc["high"] = histogram.edges.back();
    doc["observed_total"] = sum(histogram.observed);
    doc["per_library"] = nlohmann::json::array();
    for (const HistogramCheck& check : histogram.per_library) {
        nlohmann::json entry;
        entry["library"] = check.library;
        entry["expected_total"] = sum(check.expected);
        entry["skipped"] = skippedJson(check.skipped);
        if (check.merged) {
            entry["merged"] = {{"edges", check.merged->edges},
                               {"observed", check.merged->observed},
                               {"expected", check.merged->expected}};
        }
        if (check.comparison) {
            entry["chi2"] = check.comparison->chi2;
            entry["ndf"] = check.comparison->ndf;
            entry["p_value"] = check.comparison->p_value;
            entry["min_expected"] = check.comparison->min_expected;
            entry["pulls"] = check.comparison->pulls;
        }
        doc["per_library"].push_back(std::move(entry));
    }
    return doc;
}

nlohmann::json categoryObject(const std::array<double, kFirstInteractionCategories>& values) {
    nlohmann::json doc = nlohmann::json::object();
    for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
        doc[std::string{kCategoryNames.at(c)}] = values.at(c);
    }
    return doc;
}

nlohmann::json fractionsJson(const FractionResult& fractions) {
    nlohmann::json doc;
    doc["observed"] = nlohmann::json::object();
    for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
        doc["observed"][std::string{kCategoryNames.at(c)}] = fractions.observed.at(c);
    }
    doc["per_library"] = nlohmann::json::array();
    for (const FractionCheck& check : fractions.per_library) {
        nlohmann::json entry;
        entry["library"] = check.library;
        entry["expected"] = categoryObject(check.expected);
        entry["skipped"] = skippedJson(check.skipped);
        if (check.comparison) {
            entry["categories"] = nlohmann::json::array();
            for (const std::int32_t code : check.categories) {
                entry["categories"].push_back(categoryName(code));
            }
            entry["absorbing_category"] = categoryName(check.absorbing_category);
            entry["chi2"] = check.comparison->chi2;
            entry["ndf"] = check.comparison->ndf;
            entry["p_value"] = check.comparison->p_value;
            entry["fractions"] = check.comparison->fractions;
            entry["sigmas"] = check.comparison->sigmas;
            entry["pulls"] = check.comparison->pulls;
        }
        doc["per_library"].push_back(std::move(entry));
    }
    return doc;
}

nlohmann::json legendreJson(const LegendreResult& legendre) {
    nlohmann::json doc;
    doc["skipped"] = skippedJson(legendre.skipped);
    doc["per_library"] = nlohmann::json::array();
    if (legendre.moments) {
        doc["samples"] = legendre.moments->samples;
        doc["values"] = legendre.moments->values;
        doc["errors"] = legendre.moments->errors;
    }
    for (const LegendreCheck& check : legendre.per_library) {
        doc["per_library"].push_back(
            {{"library", check.library}, {"expected", check.expected}, {"pulls", check.pulls}});
    }
    return doc;
}

nlohmann::json countersJson(const ConsistencyCounters& c) {
    return {{"must_be_zero",
             {{"En_outside_source", c.En_outside_source},
              {"omega_not_beam", c.omega_not_beam},
              {"valid_np_without_scatter", c.valid_np_without_scatter},
              {"vertex_outside_box", c.vertex_outside_box},
              {"invalid_proton_energy", c.invalid_proton_energy},
              {"u_overflow", c.u_overflow},
              {"theta_p_beyond_half_pi", c.theta_p_beyond_half_pi},
              {"np_kinematics_unmatched", c.np_kinematics_unmatched},
              {"np_non_hp_outside_window", c.np_non_hp_outside_window},
              {"total", c.mustBeZeroTotal()}}},
            {"expected_rare",
             {{"np_first_without_proton", c.np_first_without_proton},
              {"np_model_ambiguous", c.np_model_ambiguous}}}};
}

nlohmann::json residualJson(const ResidualSummary& r) {
    return {{"samples", r.samples}, {"mean", r.mean}, {"rms", r.rms}, {"max_abs", r.max_abs}};
}

nlohmann::json modelMixJson(const ModelMixResult& mix) {
    nlohmann::json doc;
    doc["window_MeV"] = {kChipsElasticMinEnergy, kHpElasticMaxEnergy};
    doc["classified"] = mix.classified;
    doc["observed_non_hp"] = mix.observed_non_hp;
    doc["expected_non_hp"] = mix.expected_non_hp;
    doc["variance"] = mix.variance;
    doc["skipped"] = skippedJson(mix.skipped);
    doc["pull"] = mix.skipped.empty() ? nlohmann::json(mix.pull) : nlohmann::json(nullptr);
    doc["bins"] = mix.observed_per_bin.size();
    doc["chi2_skipped"] = skippedJson(mix.chi2_skipped);
    if (mix.chi2) {
        doc["chi2"] = *mix.chi2;
        doc["ndf"] = mix.ndf;
        doc["p_value"] = mix.p_value;
    }
    doc["Ep_residuals"] = residualJson(mix.Ep_residuals);
    return doc;
}

nlohmann::json depositCountersJson(const DepositCounters& c) {
    return {{"must_be_zero",
             {{"negative_deposit", c.negative_deposit},
              {"Evis_above_Edep", c.Evis_above_Edep},
              {"p_Evis_above_p_Edep", c.p_Evis_above_p_Edep},
              {"p_Edep_above_det_Edep", c.p_Edep_above_det_Edep},
              {"p_Edep_above_Ep_slack", c.p_Edep_above_Ep_slack},
              {"proton_without_track", c.proton_without_track},
              {"deposit_without_interaction", c.deposit_without_interaction},
              {"class_mismatch", c.class_mismatch},
              {"first_interaction_uncounted", c.first_interaction_uncounted},
              {"exit_not_on_surface", c.exit_not_on_surface},
              {"end_outside_box", c.end_outside_box},
              {"kB_zero_visible_mismatch", c.kB_zero_visible_mismatch},
              {"total", c.mustBeZeroTotal()}}},
            {"expected_rare",
             {{"p_Edep_above_Ep", c.p_Edep_above_Ep},
              {"contained_end_on_surface", c.contained_end_on_surface}}}};
}

nlohmann::json optionalJson(const std::optional<double>& value) {
    return value ? nlohmann::json(*value) : nlohmann::json(nullptr);
}

nlohmann::json residualStatsJson(const ResidualStats& s) {
    const bool any = s.samples() > 0;
    return {{"samples", s.samples()},
            {"observed_mean", any ? nlohmann::json(s.observedMean()) : nlohmann::json(nullptr)},
            {"expected_mean", any ? nlohmann::json(s.expectedMean()) : nlohmann::json(nullptr)},
            {"residual_mean", any ? nlohmann::json(s.residualMean()) : nlohmann::json(nullptr)},
            {"residual_sd", optionalJson(s.residualSd())},
            {"residual_sem", optionalJson(s.residualSem())},
            {"relative_residual", optionalJson(s.relativeResidual())},
            {"relative_sem", optionalJson(s.relativeSem())}};
}

// Empty field for an undefined value, so that a CSV reader sees a missing one.
void writeOptional(std::ostream& out, const std::optional<double>& value) {
    if (value) {
        out << *value;
    }
}

void writeMeans(std::ostream& out, const ResidualStats& s) {
    const bool any = s.samples() > 0;
    out << ',';
    writeOptional(out, any ? std::optional<double>(s.observedMean()) : std::nullopt);
    out << ',';
    writeOptional(out, any ? std::optional<double>(s.expectedMean()) : std::nullopt);
    out << ',';
    writeOptional(out, any ? std::optional<double>(s.residualMean()) : std::nullopt);
    out << ',';
    writeOptional(out, s.residualSd());
}

template <typename Check> void writeLibraryHeader(std::ostream& out, const std::vector<Check>& v) {
    for (const auto& check : v) {
        out << ",expected_" << check.library;
    }
}

} // namespace

std::string_view categoryName(std::int32_t code) {
    if (code < 0 || static_cast<std::size_t>(code) >= kFirstInteractionCategories) {
        throw std::out_of_range("no FirstInteraction category with code " + std::to_string(code));
    }
    return kCategoryNames.at(static_cast<std::size_t>(code));
}

nlohmann::json resultsJson(const CheckResults& results) {
    nlohmann::json doc;
    doc["events"] = results.events;
    doc["np_events"] = results.np_events;
    doc["np_hp_events"] = results.np_hp_events;
    doc["spectrum_u"] = histogramJson(results.spectrum_u);
    doc["depth_z"] = histogramJson(results.depth_z);
    doc["fractions"] = fractionsJson(results.fractions);
    doc["legendre"] = legendreJson(results.legendre);
    doc["Ep_residuals"] = residualJson(results.Ep_residuals);
    const OffsetSummary& offset = results.Ep_physical_offset;
    doc["Ep_physical_offset"] = residualJson(offset.absolute);
    doc["Ep_physical_offset"]["relative_mean"] = offset.relative_mean;
    doc["Ep_physical_offset"]["relative_min"] = offset.relative_min;
    doc["Ep_physical_offset"]["relative_max"] = offset.relative_max;
    doc["Ep_physical_offset"]["rounding_free_relative_min"] = offset.rounding_free_relative_min;
    doc["Ep_physical_offset"]["rounding_free_relative_max"] = offset.rounding_free_relative_max;
    doc["model_mix"] = modelMixJson(results.model_mix);
    doc["counters"] = countersJson(results.counters);
    return doc;
}

nlohmann::json depositsJson(const DepositSetup& setup, const DepositResults& results) {
    nlohmann::json doc;
    doc["skipped"] = nullptr;
    doc["setup"] = {{"half_size_mm", {setup.half_size.x, setup.half_size.y, setup.half_size.z}},
                    {"En_max_MeV", setup.En_max},
                    {"birks_kB_mm_per_MeV", setup.birks_kB / kBirksUnit},
                    {"bins", setup.bins},
                    {"deposit_sum_tolerance_MeV", kDepositSumTolerance},
                    {"nuclear_stopping_slack_MeV", kNuclearStoppingSlack},
                    {"surface_tolerance_mm", kSurfaceTolerance}};
    doc["events"] = results.events;
    doc["np_events"] = results.np_events;
    nlohmann::json classes = nlohmann::json::object();
    for (std::size_t c = 0; c < kEventClasses; ++c) {
        classes[std::string(kEventClassNames.at(c))] = results.classes.at(c);
    }
    doc["event_classes"] = classes;
    doc["events_secondary_n_hadronic"] = results.events_secondary_n_hadronic;
    doc["protons_hadronic"] = results.protons_hadronic;
    doc["protons_nuclear_recoil"] = results.protons_nuclear_recoil;
    doc["protons_exiting"] = results.protons_exiting;
    const BalanceSummary& b = results.contained_balance;
    doc["contained_balance_MeV"] = {{"samples", b.samples},   {"mean", b.mean},
                                    {"min", b.min},           {"max", b.max},
                                    {"positive", b.positive}, {"negative", b.negative}};
    doc["kB_zero_max_abs_difference_MeV"] = results.kB_zero_max_abs_difference;
    doc["Edep_overflow"] = results.Edep_overflow;
    doc["Evis_overflow"] = results.Evis_overflow;
    doc["counters"] = depositCountersJson(results.counters);
    return doc;
}

nlohmann::json depositsSkippedJson(const std::string& reason) {
    return {{"skipped", reason}};
}

nlohmann::json pstarJson(const PstarResults& results) {
    nlohmann::json bins = nlohmann::json::array();
    for (std::size_t i = 0; i < results.range_bins.size(); ++i) {
        bins.push_back({{"Ep_low_MeV", results.edges.at(i)},
                        {"Ep_high_MeV", results.edges.at(i + 1)},
                        {"range_mm", residualStatsJson(results.range_bins.at(i))},
                        {"light_MeV", residualStatsJson(results.light_bins.at(i))}});
    }
    return {{"samples", results.samples},
            {"below_table", results.below_table},
            {"range_mm", residualStatsJson(results.range)},
            {"light_MeV", residualStatsJson(results.light)},
            {"bins", bins}};
}

void writePstarComparisonCsv(std::ostream& out, const PstarResults& results) {
    setFullPrecision(out);
    out << "Ep_low,Ep_high,protons,length_mean,R_csda_mean,length_residual_mean,"
           "length_residual_sd,Evis_mean,L_mean,Evis_residual_mean,Evis_residual_sd\n";
    for (std::size_t i = 0; i < results.range_bins.size(); ++i) {
        const ResidualStats& range = results.range_bins.at(i);
        out << results.edges.at(i) << ',' << results.edges.at(i + 1) << ',' << range.samples();
        writeMeans(out, range);
        writeMeans(out, results.light_bins.at(i));
        out << '\n';
    }
}

void writeDepositSpectraCsv(std::ostream& out, const DepositResults& results) {
    setFullPrecision(out);
    out << "low,high,Edep,Evis\n";
    for (std::size_t i = 0; i < results.Edep_counts.size(); ++i) {
        out << results.edges.at(i) << ',' << results.edges.at(i + 1) << ','
            << results.Edep_counts.at(i) << ',' << results.Evis_counts.at(i) << '\n';
    }
}

void writeProtonExitsCsv(std::ostream& out, const DepositResults& results) {
    setFullPrecision(out);
    out << "Ep_low,Ep_high,protons,exiting\n";
    for (std::size_t i = 0; i < results.protons_per_Ep.size(); ++i) {
        out << results.edges.at(i) << ',' << results.edges.at(i + 1) << ','
            << results.protons_per_Ep.at(i) << ',' << results.exiting_per_Ep.at(i) << '\n';
    }
}

void writeHistogramCsv(std::ostream& out, const HistogramResult& histogram) {
    setFullPrecision(out);
    out << "low,high,observed";
    writeLibraryHeader(out, histogram.per_library);
    out << '\n';
    for (std::size_t i = 0; i < histogram.observed.size(); ++i) {
        out << histogram.edges.at(i) << ',' << histogram.edges.at(i + 1) << ','
            << histogram.observed.at(i);
        for (const HistogramCheck& check : histogram.per_library) {
            out << ',' << check.expected.at(i);
        }
        out << '\n';
    }
}

void writeFractionsCsv(std::ostream& out, const FractionResult& fractions) {
    setFullPrecision(out);
    out << "code,category,observed";
    writeLibraryHeader(out, fractions.per_library);
    out << '\n';
    for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
        out << c << ',' << kCategoryNames.at(c) << ',' << fractions.observed.at(c);
        for (const FractionCheck& check : fractions.per_library) {
            out << ',' << check.expected.at(c);
        }
        out << '\n';
    }
}

void writeLegendreCsv(std::ostream& out, const LegendreResult& legendre) {
    setFullPrecision(out);
    out << "l,value,error";
    for (const LegendreCheck& check : legendre.per_library) {
        out << ",expected_" << check.library << ",pull_" << check.library;
    }
    out << '\n';
    if (!legendre.moments) {
        return;
    }
    for (std::size_t l = 0; l < legendre.moments->values.size(); ++l) {
        out << l + 1 << ',' << legendre.moments->values.at(l) << ','
            << legendre.moments->errors.at(l);
        for (const LegendreCheck& check : legendre.per_library) {
            out << ',' << check.expected.at(l) << ',' << check.pulls.at(l);
        }
        out << '\n';
    }
}

void writeModelMixCsv(std::ostream& out, const ModelMixResult& mix) {
    setFullPrecision(out);
    out << "En_low_MeV,En_high_MeV,classified,observed_non_hp,expected_non_hp,variance\n";
    for (std::size_t i = 0; i < mix.observed_per_bin.size(); ++i) {
        out << mix.edges.at(i) << ',' << mix.edges.at(i + 1) << ',' << mix.classified_per_bin.at(i)
            << ',' << mix.observed_per_bin.at(i) << ',' << mix.expected_per_bin.at(i) << ','
            << mix.variance_per_bin.at(i) << '\n';
    }
}

void writeResidualsCsv(std::ostream& out, const std::vector<double>& mu,
                       const std::vector<double>& Ep_residuals,
                       const std::vector<double>& Ep_physical_offsets) {
    if (mu.size() != Ep_residuals.size() || mu.size() != Ep_physical_offsets.size()) {
        throw std::invalid_argument("writeResidualsCsv: the columns differ in length");
    }
    setFullPrecision(out);
    out << "mu,Ep_residual_MeV,Ep_physical_offset_MeV\n";
    for (std::size_t i = 0; i < mu.size(); ++i) {
        out << mu.at(i) << ',' << Ep_residuals.at(i) << ',' << Ep_physical_offsets.at(i) << '\n';
    }
}

void writeJsonAtomically(const std::filesystem::path& path, const nlohmann::json& doc) {
    if (std::filesystem::exists(path)) {
        throw ReportError("refusing to overwrite " + path.string());
    }
    auto temporary = path;
    temporary += kTemporarySuffix;
    {
        std::ofstream out(temporary, std::ios::out | std::ios::trunc);
        out << doc.dump(kJsonIndent) << '\n';
        out.close();
        if (!out) {
            throw ReportError("cannot write " + temporary.string());
        }
    }
    // rename is atomic within one filesystem.
    std::error_code error;
    std::filesystem::rename(temporary, path, error);
    if (error) {
        throw ReportError("cannot rename " + temporary.string() + " to " + path.string() + ": " +
                          error.message());
    }
}

} // namespace riptide::sim::checks
