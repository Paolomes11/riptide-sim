#pragma once

// Inputs of the physics checks that need Geant4 or the data root: the atom
// densities of the simulated scintillator, the H-1 target mass of the HP
// elastic model read from G4NDL, and the reference evaluations. Shared by
// riptide_sim_check (checks of a run) and riptide_sim_predict (expectations
// written before a run), so both see exactly the same medium and models.

#include <nlohmann/json.hpp>

#include <filesystem>
#include <string>
#include <vector>

#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/endf_reader.hpp"
#include "riptide/sim/hp_recoil.hpp"

namespace riptide::sim::checks {

// Atom densities of the scintillator [atoms/mm^3].
struct MediumDensities {
    double density_g_cm3 = 0.0;
    double n_H1 = 0.0;
    double n_H2 = 0.0;
    double n_C12 = 0.0;
    double n_C13 = 0.0;
};

// From the Geant4 NIST material of the simulation. Throws std::runtime_error
// if the material is missing, has an isotope without a reference evaluation,
// or lacks protium or carbon-12.
[[nodiscard]] MediumDensities scintillatorDensities();

[[nodiscard]] nlohmann::json densitiesJson(const MediumDensities& densities);

// Evaluated files under <root>/reference/<library>/, each recorded with its
// sha256 when read.
class ReferenceSet {
public:
    explicit ReferenceSet(std::filesystem::path root);

    // Throws std::runtime_error if the file is missing.
    [[nodiscard]] endf::EndfMaterial read(const std::string& library, const std::string& nuclide,
                                          const std::string& file);

    [[nodiscard]] nlohmann::json json() const;

private:
    struct File {
        std::string library;
        std::string nuclide;
        std::filesystem::path path;
    };
    std::filesystem::path root_;
    std::vector<File> files_;
};

// Header of the H-1 elastic final-state file, read as the HP model reads it.
struct HpElasticH1 {
    std::filesystem::path file;
    HpElasticHeader header;
    double target_mass = 0.0; // [MeV]
};

// Throws std::runtime_error if G4NEUTRONHPDATA is unset or the file missing.
[[nodiscard]] HpElasticH1 readHpElasticH1();

[[nodiscard]] nlohmann::json hpJson(const HpElasticH1& hp);

// JEFF-3.3, with carbon composed as the HP model composes it, then
// ENDF/B-VIII.0 with C-12 and C-13 as physics.
[[nodiscard]] std::vector<ReferenceModel> referenceModels(ReferenceSet& refs,
                                                          const MediumDensities& densities);

// Geant4 version as "major.minor.patch".
[[nodiscard]] std::string geant4Version();

} // namespace riptide::sim::checks
