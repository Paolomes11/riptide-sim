#include "riptide/sim/detector_construction.hpp"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>

#include <utility>

#include "riptide/sim/detector_materials.hpp"

namespace riptide::sim {

namespace {

constexpr const char* kWorldMaterial = "G4_Galactic";
constexpr bool kCheckOverlaps = true;

} // namespace

DetectorConstruction::DetectorConstruction(RunConfig config) : config_(std::move(config)) {}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    auto* nist = G4NistManager::Instance();
    const double world_half = worldHalfSize(config_);

    // Geant4 owns the solids, logical and physical volumes (they register
    // themselves in the geometry stores), hence the bare `new`.
    auto* world_solid = new G4Box(kWorldName, world_half, world_half, world_half);
    auto* world_logical =
        new G4LogicalVolume(world_solid, nist->FindOrBuildMaterial(kWorldMaterial), kWorldName);
    auto* world_physical = new G4PVPlacement(nullptr, G4ThreeVector(), world_logical, kWorldName,
                                             nullptr, false, 0, kCheckOverlaps);

    const auto& size = config_.geometry.scint_size;
    auto* scint_solid = new G4Box(kScintillatorName, 0.5 * size.x, 0.5 * size.y, 0.5 * size.z);
    auto* scint_logical = new G4LogicalVolume(
        scint_solid, nist->FindOrBuildMaterial(kScintillatorMaterial), kScintillatorName);
    new G4PVPlacement(nullptr, G4ThreeVector(), scint_logical, kScintillatorName, world_logical,
                      false, 0, kCheckOverlaps);
    return world_physical;
}

} // namespace riptide::sim
