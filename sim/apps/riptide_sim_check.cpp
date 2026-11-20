// Physics checks of one riptide_sim run against evaluated nuclear data: reads
// the run once, accumulates observed and expected quantities for JEFF-3.3 (the
// library Geant4's G4NDL is built from, so it checks the implementation) and
// ENDF/B-VIII.0 (an independent evaluation), and writes the numbers. It gives
// no verdict: the pass thresholds are fixed separately, before the runs. Runs
// with deposit columns are also compared with the NIST PSTAR range and
// stopping power of the scintillator.
//
// Geant4 is used only to build the scintillator material, so the atom
// densities are those of the simulated medium, and to read the H-1 target
// mass of the HP elastic model from G4NDL; no run manager is created.

#include <G4Version.hh>
#include <nlohmann/json.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/event_schema.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/version.hpp"
#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/check_cli_options.hpp"
#include "riptide/sim/check_inputs.hpp"
#include "riptide/sim/check_report.hpp"
#include "riptide/sim/deposit_checks.hpp"
#include "riptide/sim/event_row.hpp"
#include "riptide/sim/pstar.hpp"
#include "riptide/sim/pstar_checks.hpp"
#include "riptide/sim/run_config.hpp"
#include "riptide/sim/run_reader.hpp"
#include "riptide/sim/sha256.hpp"

namespace {

namespace fs = std::filesystem;
namespace checks = riptide::sim::checks;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;

// Generous upper bound on the output per simulated event: the only table that
// grows with the run is residuals_Ep.csv, one row of three full-precision
// doubles (about 75 bytes) per HP n-p event; model_mix.csv has a fixed number
// of rows.
constexpr std::uint64_t kOutputBytesPerEvent = 96;

void printVersion() {
    std::cout << "riptide_sim_check " << riptide::kProjectVersion << " (" << riptide::kGitDescribe
              << ")\n"
              << "Geant4 " << checks::geant4Version() << ' ' << G4Version << ' ' << G4Date << '\n'
              << "event schema " << riptide::schema::kSchemaVersion << '\n'
              << "checks version " << checks::kChecksVersion << '\n';
}

// The NIST PSTAR table of the vinyltoluene-based plastic, under the data root.
constexpr const char* kPstarDirectory = "pstar";
constexpr const char* kPstarFileName = "pstar_216.txt";
// The table holds for one material only: the Geant4 medium must have its
// density, up to the rounding of the material definition.
constexpr double kPstarDensityTolerance = 1e-12;

struct PstarInput {
    fs::path file;
    checks::PstarTable table;
};

PstarInput readPstar(const fs::path& root, const checks::MediumDensities& densities) {
    const double reference = checks::kVinyltolueneDensityGramPerCm3;
    if (!(std::abs(densities.density_g_cm3 - reference) <= kPstarDensityTolerance * reference)) {
        throw std::runtime_error("scintillator density " + std::to_string(densities.density_g_cm3) +
                                 " g/cm3 is not the one of the PSTAR table");
    }
    PstarInput input;
    input.file = root / "reference" / kPstarDirectory / kPstarFileName;
    std::ifstream in(input.file);
    if (!in) {
        throw std::runtime_error("reference file not found: " + input.file.string());
    }
    input.table = checks::parsePstarTable(in, reference);
    return input;
}

template <typename Writer> void writeCsv(const fs::path& path, Writer&& write) {
    if (fs::exists(path)) {
        throw std::runtime_error("refusing to overwrite " + path.string());
    }
    std::ofstream out(path);
    write(out);
    out.close();
    if (!out) {
        throw std::runtime_error("cannot write " + path.string());
    }
}

int runChecks(const checks::CheckCliOptions& opts) {
    namespace data = riptide::data;
    const fs::path root = data::resolveDataRoot(opts.data_dir);
    const fs::path run_dir = data::simRunDir(root, opts.run_id);
    const fs::path out_dir = data::validationDir(root, opts.run_id, opts.check_id);
    if (fs::exists(out_dir)) {
        throw data::DataPathError("check directory already exists: " + out_dir.string());
    }

    const nlohmann::json record = riptide::sim::readRunRecord(run_dir);
    const auto events = record.at("events").get<std::uint64_t>();
    data::requireWritableRoot(root, data::estimateRunBytes(events, kOutputBytesPerEvent),
                              data::kFreeSpaceMargin, false);
    const riptide::sim::RunConfig config = riptide::sim::parseRunConfig(record.at("config"));

    const checks::MediumDensities densities = checks::scintillatorDensities();
    const checks::HpElasticH1 hp = checks::readHpElasticH1();
    const checks::HpDeuterium deuterium = checks::readHpDeuterium();
    checks::ReferenceSet refs(root);
    checks::CheckAccumulator acc(
        checks::makeCheckSetup(config, hp.target_mass, opts.bins_u, opts.bins_z),
        checks::referenceModels(refs, densities, deuterium));
    // Runs of schema 1 (phase 1b) have no deposit columns: their checks are
    // skipped, and the truth checks above still reproduce the phase 1b numbers.
    const auto schema_version = record.at("schema_version").get<std::int32_t>();
    std::optional<checks::DepositAccumulator> deposits;
    std::optional<PstarInput> pstar_input;
    std::optional<checks::PstarComparison> pstar;
    if (schema_version >= checks::kFirstDepositSchemaVersion) {
        deposits.emplace(checks::makeDepositSetup(config));
        pstar_input = readPstar(root, densities);
        pstar.emplace(deposits->setup(), pstar_input->table);
    }
    const std::size_t class_column = riptide::sim::columnIndex("truth_event_class");
    std::cout << "riptide_sim_check: reading " << run_dir << " (" << events << " events)\n";

    const std::size_t rows =
        riptide::sim::forEachEventRow(run_dir, record, [&](const riptide::sim::EventRow& row) {
            const riptide::sim::TruthRecord truth = riptide::sim::truthRecord(row);
            acc.add(truth);
            if (deposits) {
                const riptide::sim::DepositRecord deposit = riptide::sim::depositRecord(row);
                deposits->add(truth, deposit, std::get<std::int32_t>(row.at(class_column)));
                pstar->add(truth, deposit);
            }
        });
    if (rows != events) {
        throw std::runtime_error("run " + opts.run_id + " has " + std::to_string(rows) +
                                 " rows for " + std::to_string(events) + " events");
    }
    const checks::CheckResults results = acc.finish();
    const std::optional<checks::DepositResults> deposit_results =
        deposits ? std::optional(deposits->finish()) : std::nullopt;
    const std::optional<checks::PstarResults> pstar_results =
        pstar ? std::optional(pstar->finish()) : std::nullopt;

    fs::create_directories(out_dir);
    writeCsv(out_dir / checks::kSpectrumUFileName,
             [&](std::ostream& o) { checks::writeHistogramCsv(o, results.spectrum_u); });
    writeCsv(out_dir / checks::kDepthZFileName,
             [&](std::ostream& o) { checks::writeHistogramCsv(o, results.depth_z); });
    writeCsv(out_dir / checks::kFractionsFileName,
             [&](std::ostream& o) { checks::writeFractionsCsv(o, results.fractions); });
    writeCsv(out_dir / checks::kLegendreFileName,
             [&](std::ostream& o) { checks::writeLegendreCsv(o, results.legendre); });
    writeCsv(out_dir / checks::kModelMixFileName,
             [&](std::ostream& o) { checks::writeModelMixCsv(o, results.model_mix); });
    writeCsv(out_dir / checks::kResidualsFileName, [&](std::ostream& o) {
        checks::writeResidualsCsv(o, acc.cosines(), acc.EpResiduals(), acc.EpPhysicalOffsets());
    });

    if (deposit_results) {
        writeCsv(out_dir / checks::kDepositSpectraFileName,
                 [&](std::ostream& o) { checks::writeDepositSpectraCsv(o, *deposit_results); });
        writeCsv(out_dir / checks::kProtonExitsFileName,
                 [&](std::ostream& o) { checks::writeProtonExitsCsv(o, *deposit_results); });
        writeCsv(out_dir / checks::kPstarComparisonFileName,
                 [&](std::ostream& o) { checks::writePstarComparisonCsv(o, *pstar_results); });
    }

    nlohmann::json doc;
    doc["checks_version"] = checks::kChecksVersion;
    doc["run_id"] = opts.run_id;
    doc["check_id"] = opts.check_id;
    doc["versions"] = {{"riptide", riptide::kProjectVersion},
                       {"git", riptide::kGitDescribe},
                       {"geant4", checks::geant4Version()},
                       {"schema", riptide::schema::kSchemaVersion},
                       {"run_record", record.at("run_record_version")}};
    doc["reference_files"] = refs.json();
    doc["hp_elastic_H1"] = checks::hpJson(hp);
    doc["hp_deuterium"] = checks::hpDeuteriumJson(deuterium);
    doc["medium"] = checks::densitiesJson(densities);
    doc["beam"] = {{"L_mm", acc.setup().chord.length()},
                   {"En_min_MeV", acc.setup().En_min},
                   {"En_max_MeV", acc.setup().En_max},
                   {"bins_u", opts.bins_u},
                   {"bins_z", opts.bins_z}};
    doc["results"] = checks::resultsJson(results);
    doc["results"]["deposits"] =
        deposit_results
            ? checks::depositsJson(deposits->setup(), *deposit_results)
            : checks::depositsSkippedJson("schema version " + std::to_string(schema_version) +
                                          " has no deposit columns");
    if (pstar_results) {
        nlohmann::json p = checks::pstarJson(*pstar_results);
        p["skipped"] = nullptr;
        p["file"] = pstar_input->file.string();
        p["sha256"] = checks::sha256File(pstar_input->file);
        p["density_g_cm3"] = checks::kVinyltolueneDensityGramPerCm3;
        doc["results"]["deposits"]["pstar"] = p;
    }
    // Last, and atomically: its presence marks a complete check.
    checks::writeJsonAtomically(out_dir / checks::kChecksFileName, doc);

    std::cout << "riptide_sim_check: " << results.np_events << " n-p events ("
              << results.np_hp_events << " HP), " << results.counters.mustBeZeroTotal()
              << " consistency violations";
    if (deposit_results) {
        std::cout << ", " << deposit_results->counters.mustBeZeroTotal() << " deposit violations";
    }
    std::cout << ", output in " << out_dir << '\n';
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    checks::CheckCliOptions opts;
    try {
        opts = checks::parseCheckCommandLine(std::vector<std::string>(argv + 1, argv + argc));
    } catch (const riptide::sim::CliError& e) {
        std::cerr << "riptide_sim_check: " << e.what() << "\n\n" << checks::checkHelpText();
        return kExitUsage;
    }
    if (opts.show_help) {
        std::cout << checks::checkHelpText();
        return kExitOk;
    }
    if (opts.show_version) {
        printVersion();
        return kExitOk;
    }
    try {
        return runChecks(opts);
    } catch (const std::exception& e) {
        std::cerr << "riptide_sim_check: " << e.what() << '\n';
        return kExitRuntime;
    }
}
