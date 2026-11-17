#pragma once

#include <G4VUserPrimaryGeneratorAction.hh>

#include <memory>

#include "riptide/sim/run_config.hpp"

class G4Event;
class G4ParticleGun;

namespace riptide::sim {

// One neutron per event from a fixed point along a fixed direction, with an
// energy either fixed or uniform in [En_min, En_max].
class PrimaryGeneratorAction final : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(const SourceConfig& source);
    ~PrimaryGeneratorAction() override;
    PrimaryGeneratorAction(const PrimaryGeneratorAction&) = delete;
    PrimaryGeneratorAction& operator=(const PrimaryGeneratorAction&) = delete;
    PrimaryGeneratorAction(PrimaryGeneratorAction&&) = delete;
    PrimaryGeneratorAction& operator=(PrimaryGeneratorAction&&) = delete;

    void GeneratePrimaries(G4Event* event) override;

private:
    SourceConfig source_;
    std::unique_ptr<G4ParticleGun> gun_;
};

} // namespace riptide::sim
