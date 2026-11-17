#include "riptide/sim/primary_generator_action.hpp"

#include <G4Event.hh>
#include <G4Neutron.hh>
#include <G4ParticleGun.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

namespace riptide::sim {

PrimaryGeneratorAction::PrimaryGeneratorAction(const SourceConfig& source)
    : source_(source), gun_(std::make_unique<G4ParticleGun>(1)) {
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
    if (source_.mode == SourceMode::Flat) {
        gun_->SetParticleEnergy(source_.En_min +
                                (source_.En_max - source_.En_min) * G4UniformRand());
    }
    gun_->GeneratePrimaryVertex(event);
}

} // namespace riptide::sim
