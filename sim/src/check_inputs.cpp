#include "riptide/sim/check_inputs.hpp"

#include <G4Element.hh>
#include <G4EnvironmentUtils.hh>
#include <G4Isotope.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4ParticleHPManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Version.hh>

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/sim/detector_materials.hpp"
#include "riptide/sim/sha256.hpp"

namespace riptide::sim::checks {

namespace {

namespace fs = std::filesystem;

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

void addIsotope(MediumDensities& densities, int Z, int N, double n, const std::string& name) {
    if (Z == kHydrogenZ && N == kProtium) {
        densities.n_H1 += n;
    } else if (Z == kHydrogenZ && N == kDeuterium) {
        densities.n_H2 += n;
    } else if (Z == kCarbonZ && N == kCarbon12) {
        densities.n_C12 += n;
    } else if (Z == kCarbonZ && N == kCarbon13) {
        densities.n_C13 += n;
    } else {
        throw std::runtime_error("isotope " + name + " of " + kScintillatorMaterial +
                                 " has no reference evaluation in the checks");
    }
}

} // namespace

// Geant4's internal length unit is the mm, so GetVecNbOfAtomsPerVolume is
// already in atoms/mm^3.
MediumDensities scintillatorDensities() {
    G4NistManager* nist = G4NistManager::Instance();
    const G4Material* material = nist->FindOrBuildMaterial(kScintillatorMaterial);
    if (material == nullptr) {
        throw std::runtime_error(std::string("Geant4 has no material ") + kScintillatorMaterial);
    }
    MediumDensities densities;
    densities.density_g_cm3 = material->GetDensity() / (CLHEP::g / CLHEP::cm3);
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
            addIsotope(densities, isotope->GetZ(), isotope->GetN(),
                       element_densities[e] * abundance, isotope->GetName());
        }
    }
    if (densities.n_H1 <= 0.0 || densities.n_C12 <= 0.0) {
        throw std::runtime_error("scintillator without protium or carbon-12");
    }
    return densities;
}

nlohmann::json densitiesJson(const MediumDensities& densities) {
    return {{"material", kScintillatorMaterial},
            {"density_g_cm3", densities.density_g_cm3},
            {"n_H1_per_mm3", densities.n_H1},
            {"n_H2_per_mm3", densities.n_H2},
            {"n_C12_per_mm3", densities.n_C12},
            {"n_C13_per_mm3", densities.n_C13},
            {"H2_over_H1", densities.n_H2 / densities.n_H1}};
}

ReferenceSet::ReferenceSet(fs::path root) : root_(std::move(root)) {}

endf::EndfMaterial ReferenceSet::read(const std::string& library, const std::string& nuclide,
                                      const std::string& file) {
    const fs::path path = root_ / kReferenceSubdir / library / file;
    if (!fs::is_regular_file(path)) {
        throw std::runtime_error("reference file not found: " + path.string());
    }
    files_.push_back({library, nuclide, path});
    return endf::readEndfFile(path);
}

nlohmann::json ReferenceSet::json() const {
    nlohmann::json doc = nlohmann::json::array();
    for (const File& f : files_) {
        doc.push_back({{"library", f.library},
                       {"nuclide", f.nuclide},
                       {"path", f.path.string()},
                       {"sha256", sha256File(f.path)}});
    }
    return doc;
}

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
    hp.header = parseHpElasticHeader(in);
    hp.target_mass = hp.header.awr * kinematics::kNeutronMass;
    return hp;
}

nlohmann::json hpJson(const HpElasticH1& hp) {
    return {{"file", hp.file.string()},           {"sha256", sha256File(hp.file)},
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
std::vector<ReferenceModel> referenceModels(ReferenceSet& refs, const MediumDensities& densities) {
    const endf::EndfMaterial jeff_h1 = refs.read(kJeffLibrary, "H-1", kJeffHydrogen1);
    const endf::EndfMaterial jeff_c0 = refs.read(kJeffLibrary, "C-0", kJeffCarbonNatural);
    const endf::EndfMaterial jeff_c13 = refs.read(kJeffLibrary, "C-13", kJeffCarbon13);
    ReferenceModel jeff;
    jeff.library = kJeffLibrary;
    jeff.np_angular = endf::legendreAngular(jeff_h1);
    jeff.nuclides = {makeNuclide(Target::Hydrogen, "H-1", densities.n_H1, jeff_h1),
                     makeNuclide(Target::Carbon, "C-12 (C-0)", densities.n_C12, jeff_c0),
                     makeNuclide(Target::Carbon, "C-13", densities.n_C13, jeff_c13)};

    const endf::EndfMaterial endfb_h1 = refs.read(kEndfbLibrary, "H-1", kEndfbHydrogen1);
    const endf::EndfMaterial endfb_c12 = refs.read(kEndfbLibrary, "C-12", kEndfbCarbon12);
    const endf::EndfMaterial endfb_c13 = refs.read(kEndfbLibrary, "C-13", kEndfbCarbon13);
    ReferenceModel endfb;
    endfb.library = kEndfbLibrary;
    endfb.np_angular = endf::legendreAngular(endfb_h1);
    endfb.nuclides = {makeNuclide(Target::Hydrogen, "H-1", densities.n_H1, endfb_h1),
                      makeNuclide(Target::Carbon, "C-12", densities.n_C12, endfb_c12),
                      makeNuclide(Target::Carbon, "C-13", densities.n_C13, endfb_c13)};

    std::vector<ReferenceModel> models;
    models.push_back(std::move(jeff));
    models.push_back(std::move(endfb));
    return models;
}

// G4VERSION_NUMBER packs the version as major*100 + minor*10 + patch (1141 for 11.4.1).
std::string geant4Version() {
    std::ostringstream out;
    out << G4VERSION_NUMBER / 100 << '.' << (G4VERSION_NUMBER / 10) % 10 << '.'
        << G4VERSION_NUMBER % 10;
    return out.str();
}

} // namespace riptide::sim::checks
