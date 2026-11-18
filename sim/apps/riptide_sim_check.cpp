// Physics checks of one riptide_sim run against evaluated nuclear data: reads
// the run once, accumulates observed and expected quantities for JEFF-3.3 (the
// library Geant4's G4NDL is built from, so it checks the implementation) and
// ENDF/B-VIII.0 (an independent evaluation), and writes the numbers. It gives
// no verdict: the pass thresholds are fixed separately, before the runs.
//
// Geant4 is used only to build the scintillator material, so the atom
// densities are those of the simulated medium, and to read the H-1 target
// mass of the HP elastic model from G4NDL; no run manager is created.

#include <G4Element.hh>
#include <G4EnvironmentUtils.hh>
#include <G4Isotope.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4ParticleHPManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Version.hh>
#include <nlohmann/json.hpp>

#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/event_schema.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/version.hpp"
#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/check_cli_options.hpp"
#include "riptide/sim/check_report.hpp"
#include "riptide/sim/detector_materials.hpp"
#include "riptide/sim/endf_reader.hpp"
#include "riptide/sim/event_row.hpp"
#include "riptide/sim/hp_recoil.hpp"
#include "riptide/sim/run_config.hpp"
#include "riptide/sim/run_reader.hpp"
#include "riptide/sim/sha256.hpp"

namespace {

namespace fs = std::filesystem;
namespace checks = riptide::sim::checks;
namespace endf = riptide::sim::endf;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;

// Generous upper bound on the output per simulated event: the only table that
// grows with the run is residuals_Ep.csv, one row of three full-precision
// doubles (about 75 bytes) per n-p event.
constexpr std::uint64_t kOutputBytesPerEvent = 96;

// Evaluated files under <root>/reference/<library>/, as listed with their
// sha256 in docs/architecture/data_layout.md.
constexpr const char* kReferenceSubdir = "reference";
constexpr const char* kJeffLibrary = "jeff33";
constexpr const char* kEndfbLibrary = "endfb80";
constexpr const char* kJeffHydrogen1 = "n_0125_1-H-1.dat";
constexpr const char* kJeffCarbonNatural = "n_0600_6-C-0.dat";
constexpr const char* kJeffCarbon13 = "n_0628_6-C-13.dat";
constexpr const char* kEndfbHydrogen1 = "n-001_H_001.endf";
constexpr const char* kEndfbCarbon12 = "n-006_C_012.endf";
constexpr const char* kEndfbCarbon13 = "n-006_C_013.endf";

// G4NDL elastic final state of H-1, relative to G4NEUTRONHPDATA, without the
// ".z" that G4ParticleHPManager::GetDataStream appends, as the HP model does.
constexpr const char* kHpDataEnv = "G4NEUTRONHPDATA";
constexpr const char* kHpElasticH1 = "/Elastic/FS/1_1_Hydrogen";
constexpr const char* kHpCompressedSuffix = ".z";

constexpr int kHydrogenZ = 1;
constexpr int kCarbonZ = 6;
constexpr int kProtium = 1;
constexpr int kDeuterium = 2;
constexpr int kCarbon12 = 12;
constexpr int kCarbon13 = 13;

// G4VERSION_NUMBER packs the version as major*100 + minor*10 + patch (1141 for 11.4.1).
std::string geant4Version() {
    std::ostringstream out;
    out << G4VERSION_NUMBER / 100 << '.' << (G4VERSION_NUMBER / 10) % 10 << '.'
        << G4VERSION_NUMBER % 10;
    return out.str();
}

void printVersion() {
    std::cout << "riptide_sim_check " << riptide::kProjectVersion << " (" << riptide::kGitDescribe
              << ")\n"
              << "Geant4 " << geant4Version() << ' ' << G4Version << ' ' << G4Date << '\n'
              << "event schema " << riptide::schema::kSchemaVersion << '\n'
              << "checks version " << checks::kChecksVersion << '\n';
}

// Atom densities of the scintillator [atoms/mm^3]. Geant4's internal length
// unit is the mm, so GetVecNbOfAtomsPerVolume is already in these units.
struct MediumDensities {
    double density_g_cm3 = 0.0;
    double n_H1 = 0.0;
    double n_H2 = 0.0;
    double n_C12 = 0.0;
    double n_C13 = 0.0;
};

void addIsotope(MediumDensities& d, int Z, int N, double n, const std::string& name) {
    if (Z == kHydrogenZ && N == kProtium) {
        d.n_H1 += n;
    } else if (Z == kHydrogenZ && N == kDeuterium) {
        d.n_H2 += n;
    } else if (Z == kCarbonZ && N == kCarbon12) {
        d.n_C12 += n;
    } else if (Z == kCarbonZ && N == kCarbon13) {
        d.n_C13 += n;
    } else {
        throw std::runtime_error("isotope " + name + " of " + riptide::sim::kScintillatorMaterial +
                                 " has no reference evaluation in the checks");
    }
}

MediumDensities scintillatorDensities() {
    G4NistManager* nist = G4NistManager::Instance();
    const G4Material* material = nist->FindOrBuildMaterial(riptide::sim::kScintillatorMaterial);
    if (material == nullptr) {
        throw std::runtime_error(std::string("Geant4 has no material ") +
                                 riptide::sim::kScintillatorMaterial);
    }
    MediumDensities d;
    d.density_g_cm3 = material->GetDensity() / (CLHEP::g / CLHEP::cm3);
    const G4double* element_densities = material->GetVecNbOfAtomsPerVolume();
    for (std::size_t e = 0; e < material->GetNumberOfElements(); ++e) {
        const G4Element* element = material->GetElement(static_cast<G4int>(e));
        const G4double* abundances = element->GetRelativeAbundanceVector();
        for (std::size_t i = 0; i < element->GetNumberOfIsotopes(); ++i) {
            const double abundance = abundances[i];
            const G4Isotope* isotope = element->GetIsotope(static_cast<G4int>(i));
            // Isotopes listed with zero abundance do not interact.
            if (abundance <= 0.0) {
                continue;
            }
            addIsotope(d, isotope->GetZ(), isotope->GetN(), element_densities[e] * abundance,
                       isotope->GetName());
        }
    }
    if (d.n_H1 <= 0.0 || d.n_C12 <= 0.0) {
        throw std::runtime_error("scintillator without protium or carbon-12");
    }
    return d;
}

struct ReferenceFile {
    std::string library;
    std::string nuclide;
    fs::path path;
};

class ReferenceSet {
public:
    explicit ReferenceSet(fs::path root) : root_(std::move(root)) {}

    endf::EndfMaterial read(const std::string& library, const std::string& nuclide,
                            const std::string& file) {
        const fs::path path = root_ / kReferenceSubdir / library / file;
        if (!fs::is_regular_file(path)) {
            throw std::runtime_error("reference file not found: " + path.string());
        }
        files_.push_back({library, nuclide, path});
        return endf::readEndfFile(path);
    }

    [[nodiscard]] nlohmann::json json() const {
        nlohmann::json doc = nlohmann::json::array();
        for (const ReferenceFile& f : files_) {
            doc.push_back({{"library", f.library},
                           {"nuclide", f.nuclide},
                           {"path", f.path.string()},
                           {"sha256", checks::sha256File(f.path)}});
        }
        return doc;
    }

private:
    fs::path root_;
    std::vector<ReferenceFile> files_;
};

// Header of the H-1 elastic final-state file, read as the HP model reads it.
struct HpElasticH1 {
    fs::path file;
    checks::HpElasticHeader header;
    double target_mass = 0.0; // [MeV]
};

HpElasticH1 readHpElasticH1() {
    const char* dir = G4FindDataDir(kHpDataEnv);
    if (dir == nullptr) {
        throw std::runtime_error(std::string(kHpDataEnv) + " is not set: no G4NDL to read");
    }
    const std::string base = std::string(dir) + kHpElasticH1;
    HpElasticH1 hp;
    hp.file = base + kHpCompressedSuffix;
    // GetDataStream falls back to an uncompressed file and leaves the stream
    // empty if neither exists; require the file the checks record.
    if (!fs::is_regular_file(hp.file)) {
        throw std::runtime_error("HP elastic file not found: " + hp.file.string());
    }
    std::istringstream in;
    G4ParticleHPManager::GetInstance()->GetDataStream(base, in);
    hp.header = checks::parseHpElasticHeader(in);
    hp.target_mass = hp.header.awr * riptide::kinematics::kNeutronMass;
    return hp;
}

nlohmann::json hpJson(const HpElasticH1& hp) {
    return {{"file", hp.file.string()},           {"sha256", checks::sha256File(hp.file)},
            {"rep_flag", hp.header.rep_flag},     {"awr", hp.header.awr},
            {"frame_flag", hp.header.frame_flag}, {"target_mass_MeV", hp.target_mass}};
}

// The JEFF-3.3 expectation is composed as the HP model composes carbon, so
// that it checks the implementation: G4NDL has no C-12 file and the model
// falls back to the natural-carbon one (derived from JEFF C-0) for the C-12
// share, while C-13 keeps its own file. C-13 therefore enters twice, once
// inside C-0 and once on its own, an artefact of the model that the
// expectation must reproduce. The G4NDL C-13 inelastic is the sum of the
// JEFF partial channels, not MT3 - MT102; using MT1 of C-13 leaves a known
// residual of about 3e-4 relative on the carbon nonelastic near 14 MeV.
// ENDF/B-VIII.0 has C-12 and C-13 as physics. Deuterium has no evaluation
// here: its share n_H2/n_H1 is reported as a known systematic.
std::vector<checks::ReferenceModel> referenceModels(ReferenceSet& refs, const MediumDensities& d) {
    using checks::makeNuclide;
    using checks::Target;

    const endf::EndfMaterial jeff_h1 = refs.read(kJeffLibrary, "H-1", kJeffHydrogen1);
    const endf::EndfMaterial jeff_c0 = refs.read(kJeffLibrary, "C-0", kJeffCarbonNatural);
    const endf::EndfMaterial jeff_c13 = refs.read(kJeffLibrary, "C-13", kJeffCarbon13);
    checks::ReferenceModel jeff;
    jeff.library = kJeffLibrary;
    jeff.np_angular = endf::legendreAngular(jeff_h1);
    jeff.nuclides = {makeNuclide(Target::Hydrogen, "H-1", d.n_H1, jeff_h1),
                     makeNuclide(Target::Carbon, "C-12 (C-0)", d.n_C12, jeff_c0),
                     makeNuclide(Target::Carbon, "C-13", d.n_C13, jeff_c13)};

    const endf::EndfMaterial endfb_h1 = refs.read(kEndfbLibrary, "H-1", kEndfbHydrogen1);
    const endf::EndfMaterial endfb_c12 = refs.read(kEndfbLibrary, "C-12", kEndfbCarbon12);
    const endf::EndfMaterial endfb_c13 = refs.read(kEndfbLibrary, "C-13", kEndfbCarbon13);
    checks::ReferenceModel endfb;
    endfb.library = kEndfbLibrary;
    endfb.np_angular = endf::legendreAngular(endfb_h1);
    endfb.nuclides = {makeNuclide(Target::Hydrogen, "H-1", d.n_H1, endfb_h1),
                      makeNuclide(Target::Carbon, "C-12", d.n_C12, endfb_c12),
                      makeNuclide(Target::Carbon, "C-13", d.n_C13, endfb_c13)};

    std::vector<checks::ReferenceModel> models;
    models.push_back(std::move(jeff));
    models.push_back(std::move(endfb));
    return models;
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

nlohmann::json densitiesJson(const MediumDensities& d) {
    return {{"material", riptide::sim::kScintillatorMaterial},
            {"density_g_cm3", d.density_g_cm3},
            {"n_H1_per_mm3", d.n_H1},
            {"n_H2_per_mm3", d.n_H2},
            {"n_C12_per_mm3", d.n_C12},
            {"n_C13_per_mm3", d.n_C13},
            {"H2_over_H1", d.n_H2 / d.n_H1}};
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

    const MediumDensities densities = scintillatorDensities();
    const HpElasticH1 hp = readHpElasticH1();
    ReferenceSet refs(root);
    checks::CheckAccumulator acc(
        checks::makeCheckSetup(config, hp.target_mass, opts.bins_u, opts.bins_z),
        referenceModels(refs, densities));
    std::cout << "riptide_sim_check: reading " << run_dir << " (" << events << " events)\n";

    const std::size_t rows =
        riptide::sim::forEachEventRow(run_dir, record, [&acc](const riptide::sim::EventRow& row) {
            acc.add(riptide::sim::truthRecord(row));
        });
    if (rows != events) {
        throw std::runtime_error("run " + opts.run_id + " has " + std::to_string(rows) +
                                 " rows for " + std::to_string(events) + " events");
    }
    const checks::CheckResults results = acc.finish();

    fs::create_directories(out_dir);
    writeCsv(out_dir / checks::kSpectrumUFileName,
             [&](std::ostream& o) { checks::writeHistogramCsv(o, results.spectrum_u); });
    writeCsv(out_dir / checks::kDepthZFileName,
             [&](std::ostream& o) { checks::writeHistogramCsv(o, results.depth_z); });
    writeCsv(out_dir / checks::kFractionsFileName,
             [&](std::ostream& o) { checks::writeFractionsCsv(o, results.fractions); });
    writeCsv(out_dir / checks::kLegendreFileName,
             [&](std::ostream& o) { checks::writeLegendreCsv(o, results.legendre); });
    writeCsv(out_dir / checks::kResidualsFileName, [&](std::ostream& o) {
        checks::writeResidualsCsv(o, acc.cosines(), acc.EpResiduals(), acc.EpPhysicalOffsets());
    });

    nlohmann::json doc;
    doc["checks_version"] = checks::kChecksVersion;
    doc["run_id"] = opts.run_id;
    doc["check_id"] = opts.check_id;
    doc["versions"] = {{"riptide", riptide::kProjectVersion},
                       {"git", riptide::kGitDescribe},
                       {"geant4", geant4Version()},
                       {"schema", riptide::schema::kSchemaVersion},
                       {"run_record", record.at("run_record_version")}};
    doc["reference_files"] = refs.json();
    doc["hp_elastic_H1"] = hpJson(hp);
    doc["medium"] = densitiesJson(densities);
    doc["beam"] = {{"L_mm", acc.setup().chord.length()},
                   {"En_min_MeV", acc.setup().En_min},
                   {"En_max_MeV", acc.setup().En_max},
                   {"bins_u", opts.bins_u},
                   {"bins_z", opts.bins_z}};
    doc["results"] = checks::resultsJson(results);
    // Last, and atomically: its presence marks a complete check.
    checks::writeJsonAtomically(out_dir / checks::kChecksFileName, doc);

    std::cout << "riptide_sim_check: " << results.np_events << " n-p events, "
              << results.counters.mustBeZeroTotal() << " consistency violations, output in "
              << out_dir << '\n';
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
