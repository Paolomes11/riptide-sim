#pragma once

#include <G4VUserPrimaryGeneratorAction.hh>

#include <array>
#include <memory>

#include "riptide/sim/event_seeds.hpp"
#include "riptide/sim/run_config.hpp"

class G4Event;
class G4ParticleGun;

namespace riptide::sim {

// One neutron per event from a fixed point along a fixed direction, with an
// energy either fixed or uniform in [En_min, En_max]. Each event first reseeds
// the thread's engine from (event_id, run seed), see event_seeds.hpp.
class PrimaryGeneratorAction final : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(const SourceConfig& source, std::uint64_t run_seed);
    ~PrimaryGeneratorAction() override;
    PrimaryGeneratorAction(const PrimaryGeneratorAction&) = delete;
    PrimaryGeneratorAction& operator=(const PrimaryGeneratorAction&) = delete;
    PrimaryGeneratorAction(PrimaryGeneratorAction&&) = delete;
    PrimaryGeneratorAction& operator=(PrimaryGeneratorAction&&) = delete;

    void GeneratePrimaries(G4Event* event) override;

    // Seeds of the event being generated or transported on this thread.
    [[nodiscard]] const EventSeeds& currentSeeds() const {
        return current_seeds_;
    }

private:
    SourceConfig source_;
    std::uint64_t run_seed_;
    std::unique_ptr<G4ParticleGun> gun_;
    EventSeeds current_seeds_;
    // Zero-terminated seed list for G4Random::setTheSeeds. A member, not a
    // local: MixMaxRng::setSeeds keeps the pointer (getTheSeeds returns it).
    std::array<long, 3> engine_seeds_{};
};

} // namespace riptide::sim
