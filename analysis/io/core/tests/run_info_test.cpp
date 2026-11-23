#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <optional>
#include <string>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/io/run_info.hpp"

namespace {

using nlohmann::json;
using riptide::io::HpMassSource;
using riptide::io::IoError;
using riptide::io::runInfoFromJson;

// AWR of H-1 in the G4NDL 4.7.1 elastic final-state file, as read in phase 3.
constexpr double kAwrH1 = 0.999167;
const std::filesystem::path kRunDir = "/data/sim/run_a";

json hpBlock() {
    return {{"file", "G4NDL4.7.1/Elastic/FS/1_1_Hydrogen"},
            {"sha256", "abc123"},
            {"rep_flag", 2},
            {"awr", kAwrH1},
            {"frame_flag", 2},
            {"target_mass_MeV", kAwrH1 * riptide::kinematics::kNeutronMass}};
}

json runRecord(int version) {
    json doc{{"run_record_version", version},
             {"run_id", "run_a"},
             {"git_describe", "v0.5.0"},
             {"schema_version", riptide::schema::kSchemaVersion},
             {"seed", 20260907},
             {"events", 1000},
             {"event_files", {"events_t0.root", "events_t1.root"}}};
    if (version >= riptide::io::kFirstRunRecordWithHp) {
        doc["hp_elastic_H1"] = hpBlock();
    }
    return doc;
}

json checkDoc(const std::string& run_id, const std::string& check_id) {
    return {{"run_id", run_id}, {"check_id", check_id}, {"hp_elastic_H1", hpBlock()}};
}

} // namespace

TEST_CASE("run.json v2 gives the run's own HP mass", "[run_info]") {
    const auto info = runInfoFromJson(kRunDir, runRecord(2), std::nullopt, std::nullopt);
    CHECK(info.run_id == "run_a");
    CHECK(info.events == 1000);
    CHECK(info.seed == 20260907);
    REQUIRE(info.event_files.size() == 2);
    CHECK(info.event_files[1] == kRunDir / "events_t1.root");
    CHECK(info.hp.source == HpMassSource::RunRecord);
    CHECK_THAT(info.hp.awr, Catch::Matchers::WithinULP(kAwrH1, 0));
}

TEST_CASE("run.json v1 needs a check for the HP mass", "[run_info]") {
    CHECK_THROWS_AS(runInfoFromJson(kRunDir, runRecord(1), std::nullopt, std::nullopt), IoError);
    const auto info =
        runInfoFromJson(kRunDir, runRecord(1), std::string{"c1"}, checkDoc("run_a", "c1"));
    CHECK(info.hp.source == HpMassSource::Check);
    CHECK(info.hp.check_id == "c1");
    CHECK_THAT(info.hp.target_mass_MeV,
               Catch::Matchers::WithinULP(kAwrH1 * riptide::kinematics::kNeutronMass, 0));
}

TEST_CASE("a check of another run or with another id is refused", "[run_info]") {
    CHECK_THROWS_AS(
        runInfoFromJson(kRunDir, runRecord(1), std::string{"c1"}, checkDoc("run_b", "c1")),
        IoError);
    CHECK_THROWS_AS(
        runInfoFromJson(kRunDir, runRecord(1), std::string{"c1"}, checkDoc("run_a", "c2")),
        IoError);
}

TEST_CASE("a check disagreeing with run.json v2 is refused", "[run_info]") {
    auto check = checkDoc("run_a", "c1");
    CHECK_NOTHROW(runInfoFromJson(kRunDir, runRecord(2), std::string{"c1"}, check));
    check["hp_elastic_H1"]["sha256"] = "def456";
    CHECK_THROWS_AS(runInfoFromJson(kRunDir, runRecord(2), std::string{"c1"}, check), IoError);
}

TEST_CASE("a target mass inconsistent with the AWR is refused", "[run_info]") {
    auto block = hpBlock();
    const double relative_shift = 1e-9;
    block["target_mass_MeV"] = block["target_mass_MeV"].get<double>() * (1.0 + relative_shift);
    CHECK_THROWS_AS(riptide::io::hpMassFromBlock(block, HpMassSource::RunRecord, {}), IoError);
    block = hpBlock();
    block["awr"] = -kAwrH1;
    CHECK_THROWS_AS(riptide::io::hpMassFromBlock(block, HpMassSource::RunRecord, {}), IoError);
    block = hpBlock();
    block.erase("target_mass_MeV");
    CHECK_THROWS_AS(riptide::io::hpMassFromBlock(block, HpMassSource::RunRecord, {}), IoError);
}

TEST_CASE("unknown versions and mismatched ids are refused", "[run_info]") {
    CHECK_THROWS_AS(runInfoFromJson(kRunDir, runRecord(0), std::nullopt, std::nullopt), IoError);
    CHECK_THROWS_AS(runInfoFromJson(kRunDir, runRecord(3), std::nullopt, std::nullopt), IoError);
    auto record = runRecord(2);
    record["schema_version"] = riptide::schema::kSchemaVersion + 1;
    CHECK_THROWS_AS(runInfoFromJson(kRunDir, record, std::nullopt, std::nullopt), IoError);
    record = runRecord(2);
    record["run_id"] = "run_b";
    CHECK_THROWS_AS(runInfoFromJson(kRunDir, record, std::nullopt, std::nullopt), IoError);
}

TEST_CASE("event file names cannot leave the run directory", "[run_info]") {
    auto record = runRecord(2);
    record["event_files"] = {"../other/events_t0.root"};
    CHECK_THROWS_AS(runInfoFromJson(kRunDir, record, std::nullopt, std::nullopt), IoError);
    record["event_files"] = json::array();
    CHECK_THROWS_AS(runInfoFromJson(kRunDir, record, std::nullopt, std::nullopt), IoError);
}

namespace {

using riptide::io::sourceSpectrumFromRunRecord;
using riptide::io::SpectrumMode;

json withSource(json source) {
    return {{"config", {{"source", std::move(source)}}}};
}

} // namespace

TEST_CASE("config.source gives the spectrum of the run", "[run_info][source]") {
    const auto mono = sourceSpectrumFromRunRecord(withSource({{"mode", "mono"}, {"En_MeV", 14.0}}));
    CHECK(mono.mode == SpectrumMode::Mono);
    CHECK_THAT(mono.En_MeV, Catch::Matchers::WithinULP(14.0, 0));

    const auto flat = sourceSpectrumFromRunRecord(
        withSource({{"mode", "flat"}, {"En_min_MeV", 10.0}, {"En_max_MeV", 19.5}}));
    CHECK(flat.mode == SpectrumMode::Flat);
    CHECK_THAT(flat.En_min_MeV, Catch::Matchers::WithinULP(10.0, 0));
    CHECK_THAT(flat.En_max_MeV, Catch::Matchers::WithinULP(19.5, 0));
    CHECK(riptide::io::spectrumModeName(flat.mode) == "flat");
}

TEST_CASE("config.source refuses what it does not know", "[run_info][source]") {
    CHECK_THROWS_AS(sourceSpectrumFromRunRecord(json::object()), IoError);
    CHECK_THROWS_AS(sourceSpectrumFromRunRecord({{"config", json::object()}}), IoError);
    CHECK_THROWS_AS(sourceSpectrumFromRunRecord(withSource({{"mode", "gauss"}, {"En_MeV", 14.0}})),
                    IoError);
    CHECK_THROWS_AS(sourceSpectrumFromRunRecord(withSource({{"mode", "mono"}})), IoError);
    CHECK_THROWS_AS(sourceSpectrumFromRunRecord(withSource({{"mode", "mono"}, {"En_MeV", -1.0}})),
                    IoError);
    CHECK_THROWS_AS(sourceSpectrumFromRunRecord(withSource({{"mode", "mono"}, {"En_MeV", "14"}})),
                    IoError);
    CHECK_THROWS_AS(
        sourceSpectrumFromRunRecord(withSource({{"mode", "flat"}, {"En_min_MeV", 10.0}})), IoError);
    CHECK_THROWS_AS(sourceSpectrumFromRunRecord(
                        withSource({{"mode", "flat"}, {"En_min_MeV", 20.0}, {"En_max_MeV", 10.0}})),
                    IoError);
    CHECK_THROWS_AS(sourceSpectrumFromRunRecord(
                        withSource({{"mode", "flat"}, {"En_min_MeV", 10.0}, {"En_max_MeV", 10.0}})),
                    IoError);
}

TEST_CASE("the source of the shipped configs is understood", "[run_info][source]") {
    // run.json copies the config file under "config", so a shipped config
    // nested the same way is what a real run gives the reader.
    const std::filesystem::path config_dir = RIPTIDE_CONFIG_DIR;
    const auto mono = sourceSpectrumFromRunRecord(
        {{"config", riptide::io::readJsonFile(config_dir / "3_ff_mono_14MeV_x.json")}});
    CHECK(mono.mode == SpectrumMode::Mono);
    CHECK_THAT(mono.En_MeV, Catch::Matchers::WithinULP(14.0, 0));
    const auto flat = sourceSpectrumFromRunRecord(
        {{"config",
          riptide::io::readJsonFile(config_dir / "3_ff_flat_10_19p5MeV_isotropic.json")}});
    CHECK(flat.mode == SpectrumMode::Flat);
    CHECK_THAT(flat.En_max_MeV, Catch::Matchers::WithinULP(19.5, 0));
}
