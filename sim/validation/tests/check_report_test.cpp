#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "riptide/sim/check_report.hpp"

namespace {

namespace fs = std::filesystem;
using Catch::Matchers::WithinULP;
using namespace riptide::sim::checks;

constexpr double kMuA = 0.25;
constexpr double kMuB = -0.125;
constexpr double kResidualA = 1.0e-3;
constexpr double kResidualB = -2.5e-4;

HistogramResult twoBinHistogram() {
    HistogramResult h;
    h.edges = {0.0, 0.5, 1.0};
    h.observed = {3.0, 7.0};
    HistogramCheck compared;
    compared.library = "jeff33";
    compared.expected = {4.0, 6.0};
    compared.merged = MergedBins{{3.0, 7.0}, {4.0, 6.0}, {0.0, 0.5, 1.0}};
    compared.comparison = HistogramComparison{};
    compared.comparison->chi2 = 0.4;
    compared.comparison->ndf = 1;
    compared.comparison->p_value = 0.5;
    compared.comparison->min_expected = 4.0;
    compared.comparison->pulls = {-0.5, 0.4};
    HistogramCheck skipped;
    skipped.library = "endfb80";
    skipped.expected = {1.0, 2.0};
    skipped.skipped = "too few expected events";
    h.per_library = {compared, skipped};
    return h;
}

CheckResults sampleResults() {
    CheckResults r;
    r.events = 100;
    r.np_events = 10;
    r.spectrum_u = twoBinHistogram();
    r.depth_z = twoBinHistogram();
    r.fractions.observed = {0, 10, 5, 3, 82};
    FractionCheck f;
    f.library = "jeff33";
    f.expected = {0.0, 11.0, 4.0, 3.0, 82.0};
    f.categories = {1, 2, 3, 4};
    f.absorbing_category = 4;
    f.comparison = FractionComparison{};
    f.comparison->total = 100;
    r.fractions.per_library = {f};
    r.legendre.moments = LegendreMoments{10, {0.1, 0.2}, {0.05, 0.06}};
    r.legendre.per_library = {LegendreCheck{"jeff33", {0.0, 0.2}, {2.0, 0.0}}};
    r.Ep_residuals = {2, 3.75e-4, 7.3e-4, 1.0e-3};
    r.counters.omega_not_beam = 2;
    r.counters.u_overflow = 1;
    return r;
}

std::vector<std::string> lines(const std::string& text) {
    std::vector<std::string> out;
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        out.push_back(line);
    }
    return out;
}

fs::path freshTempDir(const std::string& name) {
    const fs::path dir = fs::temp_directory_path() / name;
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

} // namespace

TEST_CASE("categoryName follows the FirstInteraction codes", "[check_report]") {
    CHECK(categoryName(0) == "unknown");
    CHECK(categoryName(1) == "np_elastic");
    CHECK(categoryName(2) == "nc_elastic");
    CHECK(categoryName(3) == "nc_inelastic");
    CHECK(categoryName(4) == "other");
    CHECK_THROWS_AS(categoryName(-1), std::out_of_range);
    CHECK_THROWS_AS(categoryName(5), std::out_of_range);
}

TEST_CASE("resultsJson reports compared and skipped checks", "[check_report]") {
    const nlohmann::json doc = resultsJson(sampleResults());

    CHECK(doc["events"] == 100);
    CHECK(doc["np_events"] == 10);

    const auto& u = doc["spectrum_u"];
    CHECK(u["bins"] == 2);
    CHECK_THAT(u["high"].get<double>(), WithinULP(1.0, 0));
    CHECK_THAT(u["observed_total"].get<double>(), WithinULP(10.0, 0));
    const auto& jeff = u["per_library"][0];
    CHECK(jeff["library"] == "jeff33");
    CHECK(jeff["skipped"].is_null());
    CHECK(jeff["ndf"] == 1);
    CHECK(jeff["merged"]["edges"].size() == 3);
    CHECK(jeff["pulls"].size() == 2);
    const auto& endf = u["per_library"][1];
    CHECK(endf["skipped"] == "too few expected events");
    CHECK_FALSE(endf.contains("chi2"));
    CHECK_FALSE(endf.contains("merged"));
    CHECK_THAT(endf["expected_total"].get<double>(), WithinULP(3.0, 0));

    const auto& fr = doc["fractions"];
    CHECK(fr["observed"]["other"] == 82);
    CHECK(fr["per_library"][0]["absorbing_category"] == "other");
    CHECK(fr["per_library"][0]["categories"][0] == "np_elastic");
    CHECK_THAT(fr["per_library"][0]["expected"]["nc_elastic"].get<double>(), WithinULP(4.0, 0));

    const auto& leg = doc["legendre"];
    CHECK(leg["samples"] == 10);
    CHECK(leg["values"].size() == 2);
    CHECK(leg["per_library"][0]["pulls"].size() == 2);
    CHECK(leg["skipped"].is_null());

    CHECK(doc["Ep_residuals"]["samples"] == 2);
    CHECK(doc["counters"]["must_be_zero"]["omega_not_beam"] == 2);
    CHECK(doc["counters"]["must_be_zero"]["total"] == 2);
    CHECK(doc["counters"]["expected_rare"]["u_overflow"] == 1);
}

TEST_CASE("resultsJson of a skipped Legendre check has no moments", "[check_report]") {
    CheckResults r = sampleResults();
    r.legendre = LegendreResult{};
    r.legendre.skipped = "no n-p events";
    const nlohmann::json doc = resultsJson(r);
    CHECK(doc["legendre"]["skipped"] == "no n-p events");
    CHECK_FALSE(doc["legendre"].contains("values"));
    CHECK(doc["legendre"]["per_library"].empty());
}

TEST_CASE("histogram CSV has one row per original bin", "[check_report]") {
    std::ostringstream out;
    writeHistogramCsv(out, twoBinHistogram());
    const auto rows = lines(out.str());
    REQUIRE(rows.size() == 3);
    CHECK(rows[0] == "low,high,observed,expected_jeff33,expected_endfb80");
    CHECK(rows[1] == "0,0.5,3,4,1");
    CHECK(rows[2] == "0.5,1,7,6,2");
}

TEST_CASE("fractions CSV lists every category", "[check_report]") {
    std::ostringstream out;
    writeFractionsCsv(out, sampleResults().fractions);
    const auto rows = lines(out.str());
    REQUIRE(rows.size() == 6);
    CHECK(rows[0] == "code,category,observed,expected_jeff33");
    CHECK(rows[1] == "0,unknown,0,0");
    CHECK(rows[5] == "4,other,82,82");
}

TEST_CASE("Legendre CSV starts at l = 1, header only when skipped", "[check_report]") {
    std::ostringstream out;
    writeLegendreCsv(out, sampleResults().legendre);
    const auto rows = lines(out.str());
    REQUIRE(rows.size() == 3);
    CHECK(rows[0] == "l,value,error,expected_jeff33,pull_jeff33");
    CHECK(rows[1].rfind("1,", 0) == 0);
    CHECK(rows[2].rfind("2,", 0) == 0);

    std::ostringstream skipped;
    writeLegendreCsv(skipped, LegendreResult{});
    CHECK(lines(skipped.str()) == std::vector<std::string>{"l,value,error"});
}

TEST_CASE("residuals CSV keeps every double exactly", "[check_report]") {
    const double third = 1.0 / 3.0;
    std::ostringstream out;
    writeResidualsCsv(out, {kMuA, third}, {kResidualA, kResidualB});
    const auto rows = lines(out.str());
    REQUIRE(rows.size() == 3);
    CHECK(rows[0] == "mu,Ep_residual_MeV");
    const std::string mu_text = rows[2].substr(0, rows[2].find(','));
    CHECK_THAT(std::stod(mu_text), WithinULP(third, 0));

    CHECK_THROWS_AS(writeResidualsCsv(out, {kMuA, kMuB}, {kResidualA}), std::invalid_argument);
}

TEST_CASE("writeJsonAtomically writes once and refuses to overwrite", "[check_report]") {
    const fs::path dir = freshTempDir("riptide_check_report_test");
    const fs::path path = dir / kChecksFileName;
    const nlohmann::json doc = {{"answer", 42}};

    writeJsonAtomically(path, doc);
    std::ifstream in(path);
    CHECK(nlohmann::json::parse(in) == doc);
    fs::path temporary = path;
    temporary += ".tmp";
    CHECK_FALSE(fs::exists(temporary));

    CHECK_THROWS_AS(writeJsonAtomically(path, doc), ReportError);
    CHECK_THROWS_AS(writeJsonAtomically(dir / "missing" / kChecksFileName, doc), ReportError);
    fs::remove_all(dir);
}
