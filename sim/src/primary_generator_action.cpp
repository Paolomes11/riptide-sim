#include "riptide/sim/primary_generator_action.hpp"

#include <G4CrossSectionDataStore.hh>
#include <G4DynamicParticle.hh>
#include <G4Event.hh>
#include <G4HadronicProcess.hh>
#include <G4LogicalVolume.hh>
#include <G4Navigator.hh>
#include <G4Neutron.hh>
#include <G4ParticleGun.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>
#include <G4ThreeVector.hh>
#include <G4TransportationManager.hh>
#include <G4VPhysicalVolume.hh>
#include <Randomize.hh>

namespace riptide::sim {

namespace {

// The dummy evaluation below uses an energy no primary can have: every gun
// energy is >= En_min > 0.
constexpr double kCacheResetEnergyFraction = 0.5;

// G4CrossSectionDataStore::GetCrossSection reuses the last value when the
// (energy, material, particle) key is unchanged. The HP neutron data compute
// that value by a Monte Carlo Doppler integration that draws random numbers,
// so a cache hit skips draws. The key left by the previous event on the same
// thread would then decide how many numbers the first steps consume, and the
// event would depend on the thread's history, not only on its seeds. A dummy
// evaluation at an energy outside the source range moves every key away, so
// each event starts with a miss. It runs before the reseed: its own draws do
// not belong to the event.
void resetNeutronCrossSectionCaches(double dummy_energy) {
    const auto* world = G4TransportationManager::GetTransportationManager()
                            ->GetNavigatorForTracking()
                            ->GetWorldVolume();
    const auto* material = world->GetLogicalVolume()->GetMaterial();
    const G4DynamicParticle neutron(G4Neutron::Definition(), G4ThreeVector(0.0, 0.0, 1.0),
                                    dummy_energy);
    const auto* manager = G4Neutron::Definition()->GetProcessManager();
    const auto* processes = manager->GetProcessList();
    for (G4int i = 0; i < manager->GetProcessListLength(); ++i) {
        auto* hadronic = dynamic_cast<G4HadronicProcess*>((*processes)[i]);
        if (hadronic != nullptr) {
            hadronic->GetCrossSectionDataStore()->ComputeCrossSection(&neutron, material);
        }
    }
}

} // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction(const SourceConfig& source,
                                               const GeometryConfig& geometry,
                                               std::uint64_t run_seed, TruthRecorder& recorder)
    : source_(source), run_seed_(run_seed), recorder_(&recorder),
      gun_(std::make_unique<G4ParticleGun>(1)) {
    if (source_.beam == BeamMode::FarField) {
        disk_ = farFieldDisk(source_, scintHalfSize(geometry));
    }
    gun_->SetParticleDefinition(G4Neutron::Definition());
    const auto& d = source_.direction;
    const auto& p = source_.position;
    // Project units equal Geant4's internal units (checked in sim/tests), so
    // the values are passed through unchanged.
    gun_->SetParticleMomentumDirection(G4ThreeVector(d.x, d.y, d.z));
    gun_->SetParticlePosition(G4ThreeVector(p.x, p.y, p.z));
    gun_->SetParticleEnergy(source_.En_min);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    resetNeutronCrossSectionCaches(kCacheResetEnergyFraction * source_.En_min);
    // G4WorkerRunManager::GenerateEvent reseeds the engine with its own seeds
    // and calls this method next, with no draw in between; reseeding here, before
    // the first draw of the event, replaces those seeds with ours.
    current_seeds_ = eventSeeds(run_seed_, event->GetEventID());
    engine_seeds_ = {static_cast<long>(current_seeds_.seed_0),
                     static_cast<long>(current_seeds_.seed_1), 0};
    G4Random::setTheSeeds(engine_seeds_.data());

    if (source_.mode == SourceMode::Flat) {
        gun_->SetParticleEnergy(source_.En_min +
                                (source_.En_max - source_.En_min) * G4UniformRand());
    }
    Vec3 direction = source_.direction;
    Vec3 start = source_.position;
    if (source_.beam == BeamMode::FarField) {
        const auto ray = sampleFarField(source_, disk_, [] { return G4UniformRand(); });
        direction = ray.direction;
        start = ray.position;
        gun_->SetParticleMomentumDirection(G4ThreeVector(direction.x, direction.y, direction.z));
        gun_->SetParticlePosition(G4ThreeVector(ray.position.x, ray.position.y, ray.position.z));
    }
    const EventMeta meta{event->GetEventID(), current_seeds_.seed_0, current_seeds_.seed_1};
    recorder_->beginEvent(meta, gun_->GetParticleEnergy(), direction, start);
    gun_->GeneratePrimaryVertex(event);
}

} // namespace riptide::sim
