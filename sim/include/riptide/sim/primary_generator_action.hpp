#pragma once

#include <G4VUserPrimaryGeneratorAction.hh>

#include <array>
#include <memory>

#include "riptide/sim/event_seeds.hpp"
#include "riptide/sim/run_config.hpp"
#include "riptide/sim/source_sampling.hpp"
#include "riptide/sim/truth_record.hpp"

class G4Event;
class G4ParticleGun;

namespace riptide::sim {

// One neutron per event, with an energy either fixed or uniform in
// [En_min, En_max]. A pencil beam starts from a fixed point along a fixed
// direction; a far-field source draws the direction and the start point of
// every event (D27, source_sampling.hpp) after the energy, so the pencil
// events keep the draws of the earlier phases. Each event first clears
// the neutron cross-section caches left by the previous event, then reseeds the
// thread's engine from (event_id, run seed), see event_seeds.hpp, and then
// opens the event's truth record: Geant4 calls GeneratePrimaries before
// BeginOfEventAction, so this is the earliest point of the event.
class PrimaryGeneratorAction final : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(const SourceConfig& source, const GeometryConfig& geometry,
                           std::uint64_t run_seed, TruthRecorder& recorder);
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
    FarFieldDisk disk_; // used by far-field sources only
    std::uint64_t run_seed_;
    TruthRecorder* recorder_; // owned by the thread's RunAction
    std::unique_ptr<G4ParticleGun> gun_;
    EventSeeds current_seeds_;
    // Zero-terminated seed list for G4Random::setTheSeeds. A member, not a
    // local: MixMaxRng::setSeeds keeps the pointer (getTheSeeds returns it).
    std::array<long, 3> engine_seeds_{};
};

} // namespace riptide::sim
