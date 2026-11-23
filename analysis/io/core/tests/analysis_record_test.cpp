#include <catch2/catch_test_macros.hpp>

#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>

#include "riptide/io/analysis_record.hpp"

namespace {

using riptide::io::AnalysisSettings;
using riptide::io::CutFlow;
using riptide::io::IoError;
using riptide::io::Selection;

riptide::io::RunInfo sampleRun() {
    riptide::io::RunInfo run;
    run.run_id = "run_a";
    run.run_record_version = 1;
    run.schema_version = 4;
    run.events = 10;
    run.git_describe = "v0.5.0";
    run.hp.awr = 0.999167;
    run.hp.target_mass_MeV = 938.7827020551201;
    run.hp.source = riptide::io::HpMassSource::Check;
    run.hp.check_id = "c1";
    return run;
}

AnalysisSettings sampleSettings() {
    AnalysisSettings settings;
    settings.analysis_id = "a1";
    settings.seed = 20260907;
    settings.resolution = {0.10, 0.08};
    settings.selection = Selection::parse("single_np,contained");
    return settings;
}

CutFlow sampleFlow() {
    CutFlow flow(2);
    flow.record(std::nullopt);
    flow.record(std::optional<std::size_t>{1});
    flow.record(std::optional<std::size_t>{0});
    flow.record(std::nullopt);
    return flow;
}

std::string readAll(const std::filesystem::path& path) {
    std::ifstream in(path);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

} // namespace

TEST_CASE("analysis.json records the inputs and the cut flow", "[analysis_record]") {
    const auto doc = riptide::io::analysisRecordJson(sampleRun(), sampleSettings(), sampleFlow(),
                                                     nlohmann::json::object());
    CHECK(doc.at("analysis_record_version") == riptide::io::kAnalysisRecordVersion);
    CHECK(doc.at("input").at("run_id") == "run_a");
    CHECK(doc.at("hp_target_mass").at("source") == "checks.json");
    CHECK(doc.at("hp_target_mass").at("check_id") == "c1");
    CHECK(doc.at("response").at("seed") == 20260907);
    CHECK(doc.at("selection").at("text") == "single_np,contained");
    CHECK(doc.at("selection").at("events_read") == 4);
    CHECK(doc.at("selection").at("events_selected") == 2);
    CHECK(doc.at("selection").at("cuts").at(0).at("passed") == 3);
    // Same inputs, same text: the record is reproducible.
    CHECK(doc.dump() == riptide::io::analysisRecordJson(sampleRun(), sampleSettings(), sampleFlow(),
                                                        nlohmann::json::object())
                            .dump());
}

TEST_CASE("the cut-flow CSV has one line per cut", "[analysis_record]") {
    const auto csv = riptide::io::cutFlowCsv(sampleSettings().selection, sampleFlow());
    CHECK(csv == "index,cut,threshold_MeV,passed,fraction\n"
                 "0,all,,4,1\n"
                 "1,single_np,,3,0.75\n"
                 "2,contained,,2,0.5\n");
    CHECK_THROWS_AS(riptide::io::cutFlowCsv(sampleSettings().selection, CutFlow(1)), IoError);
}

TEST_CASE("the cuts block lists each cut with its passing events", "[analysis_record]") {
    const auto cuts = riptide::io::cutsJson(sampleSettings().selection, sampleFlow());
    REQUIRE(cuts.size() == 2);
    CHECK(cuts.at(0).at("name") == "single_np");
    CHECK(cuts.at(0).at("passed") == 3);
    CHECK(cuts.at(1).at("name") == "contained");
    CHECK(cuts.at(1).at("passed") == 2);
    CHECK_FALSE(cuts.at(1).contains("threshold_MeV"));
    CHECK_THROWS_AS(riptide::io::cutsJson(sampleSettings().selection, CutFlow(1)), IoError);
}

TEST_CASE("files are written once and never overwritten", "[analysis_record]") {
    const auto dir = std::filesystem::temp_directory_path() /
                     ("riptide_io_core_test_" + std::to_string(::getpid()));
    std::filesystem::create_directories(dir);
    const auto target = dir / "analysis.json";
    riptide::io::writeFileAtomically(target, "first\n");
    CHECK(readAll(target) == "first\n");
    CHECK_FALSE(std::filesystem::exists(dir / "analysis.json.tmp"));
    CHECK_THROWS_AS(riptide::io::writeFileAtomically(target, "second\n"), IoError);
    CHECK(readAll(target) == "first\n");
    std::filesystem::remove_all(dir);
}
