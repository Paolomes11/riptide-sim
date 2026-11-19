// Expected outcome of the physics checks of a planned run, written before the
// run exists: from the run configuration, the evaluated files and the medium
// of the simulation, it prints as JSON the expected first-interaction
// fractions and counts, the expected Legendre moments and u spectrum, their
// statistical errors at the planned event count, the separation between the
// reference libraries, the range of the HP proton-energy offset and the
// acceptance plan of the statistical tests. Nothing is simulated or written to
// disk; the output is meant to be committed with the predictions.
//
// Geant4 is used, as in riptide_sim_check, only for the scintillator material
// and the H-1 target mass of the HP elastic model.

#include <G4Version.hh>
#include <nlohmann/json.hpp>

#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/version.hpp"
#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/check_cli_options.hpp"
#include "riptide/sim/check_inputs.hpp"
#include "riptide/sim/check_report.hpp"
#include "riptide/sim/run_config.hpp"
#include "riptide/sim/run_prediction.hpp"

namespace {

namespace fs = std::filesystem;
namespace checks = riptide::sim::checks;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;

// Family-wise false-alarm probability of the statistical tests of one run:
// with three runs, a correct simulation fails some test with probability
// about 3 %, small enough that a failure calls for a cause, large enough
// that the per-test thresholds stay near the usual 3-4 sigma.
constexpr double kAlphaFamily = 0.01;

// The quadrature error is estimated against twice the nodes.
constexpr std::size_t kNodeRefinement = 2;

// Output of nlohmann::json::dump, for a diff-friendly committed file.
constexpr int kJsonIndent = 2;

void printVersion() {
    std::cout << "riptide_sim_predict " << riptide::kProjectVersion << " (" << riptide::kGitDescribe
              << ")\n"
              << "Geant4 " << checks::geant4Version() << ' ' << G4Version << ' ' << G4Date << '\n'
              << "checks version " << checks::kChecksVersion << '\n';
}

int predict(const checks::PredictCliOptions& opts) {
    const fs::path root = riptide::data::resolveDataRoot(opts.data_dir);
    const riptide::sim::RunConfig config = riptide::sim::loadRunConfig(opts.config_path);

    const checks::MediumDensities densities = checks::scintillatorDensities();
    const checks::HpElasticH1 hp = checks::readHpElasticH1();
    const checks::HpDeuterium deuterium = checks::readHpDeuterium();
    checks::ReferenceSet refs(root);
    const std::vector<checks::ReferenceModel> models =
        checks::referenceModels(refs, densities, deuterium);
    const checks::CheckSetup setup =
        checks::makeCheckSetup(config, hp.target_mass, opts.bins_u, opts.bins_z);

    const checks::RunPrediction run =
        checks::predictRun(setup, models, opts.events, opts.energy_nodes);
    const checks::RunPrediction refined =
        checks::predictRun(setup, models, opts.events, kNodeRefinement * opts.energy_nodes);
    // The coarse rule against the refined one, in units of the expected errors
    // of the first library: what the quadrature adds to every comparison.
    const checks::PredictionShift quadrature =
        checks::predictionShift(refined.libraries.front(), run.libraries.front(),
                                checks::expectedErrors(refined, 0), opts.events);
    const checks::TestPlan plan = checks::testPlan(run, kAlphaFamily);

    nlohmann::json doc;
    doc["checks_version"] = checks::kChecksVersion;
    doc["config"] = opts.config_path;
    doc["versions"] = {{"riptide", riptide::kProjectVersion},
                       {"git", riptide::kGitDescribe},
                       {"geant4", checks::geant4Version()}};
    doc["reference_files"] = refs.json();
    doc["hp_elastic_H1"] = checks::hpJson(hp);
    doc["hp_deuterium"] = checks::hpDeuteriumJson(deuterium);
    doc["medium"] = checks::densitiesJson(densities);
    doc["bins_u"] = opts.bins_u;
    doc["bins_z"] = opts.bins_z;
    doc["prediction"] = checks::predictionJson(run, quadrature, plan);
    std::cout << doc.dump(kJsonIndent) << '\n';
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    checks::PredictCliOptions opts;
    try {
        opts = checks::parsePredictCommandLine(std::vector<std::string>(argv + 1, argv + argc));
    } catch (const riptide::sim::CliError& e) {
        std::cerr << "riptide_sim_predict: " << e.what() << "\n\n" << checks::predictHelpText();
        return kExitUsage;
    }
    if (opts.show_help) {
        std::cout << checks::predictHelpText();
        return kExitOk;
    }
    if (opts.show_version) {
        printVersion();
        return kExitOk;
    }
    try {
        return predict(opts);
    } catch (const std::exception& e) {
        std::cerr << "riptide_sim_predict: " << e.what() << '\n';
        return kExitRuntime;
    }
}
