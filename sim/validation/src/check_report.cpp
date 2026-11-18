#include "riptide/sim/check_report.hpp"

#include <array>
#include <cstddef>
#include <fstream>
#include <ios>
#include <limits>
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
              {"total", c.mustBeZeroTotal()}}},
            {"expected_rare",
             {{"np_first_without_proton", c.np_first_without_proton},
              {"u_overflow", c.u_overflow},
              {"theta_p_beyond_half_pi", c.theta_p_beyond_half_pi}}}};
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
    doc["spectrum_u"] = histogramJson(results.spectrum_u);
    doc["depth_z"] = histogramJson(results.depth_z);
    doc["fractions"] = fractionsJson(results.fractions);
    doc["legendre"] = legendreJson(results.legendre);
    doc["Ep_residuals"] = {{"samples", results.Ep_residuals.samples},
                           {"mean", results.Ep_residuals.mean},
                           {"rms", results.Ep_residuals.rms},
                           {"max_abs", results.Ep_residuals.max_abs}};
    doc["counters"] = countersJson(results.counters);
    return doc;
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

void writeResidualsCsv(std::ostream& out, const std::vector<double>& mu,
                       const std::vector<double>& Ep_residuals) {
    if (mu.size() != Ep_residuals.size()) {
        throw std::invalid_argument("writeResidualsCsv: mu and Ep residuals differ in length");
    }
    setFullPrecision(out);
    out << "mu,Ep_residual_MeV\n";
    for (std::size_t i = 0; i < mu.size(); ++i) {
        out << mu.at(i) << ',' << Ep_residuals.at(i) << '\n';
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
