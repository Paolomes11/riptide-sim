#pragma once

// Energy deposit and visible (Birks-quenched) energy of one event, plus the
// track of the recoil proton of the first n-p scatter. Geant4-free like
// truth_record.hpp: the Geant4 side computes edep and evis of each step in the
// scintillator and hands them over as a DepositStep.

#include <cstdint>

#include "riptide/sim/run_config.hpp"

namespace riptide::sim {

// Which track a step belongs to. Only the recoil proton of the first n-p
// scatter is followed on its own; everything else (its delta electrons
// included) only adds to the event totals.
enum class TrackRole { RecoilProton, Other };

struct DepositStep {
    TrackRole role = TrackRole::Other;
    double edep = 0.0;   // deposited energy [MeV]
    double evis = 0.0;   // visible energy after Birks quenching [MeV]
    Vec3 pre;            // pre-step point [mm]
    Vec3 post;           // post-step point [mm]
    double length = 0.0; // step length [mm]
    // The post-step point is on the scintillator boundary and the track goes
    // on outside it.
    bool leaves_volume = false;
    // The step ended in a hadronic interaction of this track.
    bool hadronic = false;
    // A non-hadronic step (Coulomb scattering) whose secondaries include a
    // nucleus or a hydrogen proton: that recoil takes energy off the track.
    bool nuclear_recoil = false;
};

// Track of the recoil proton inside the scintillator.
struct ProtonTrack {
    bool seen = false;           // at least one step recorded; the fields below are 0 otherwise
    Vec3 start;                  // pre-step point of its first step [mm]
    Vec3 end;                    // post-step point of its last step [mm]
    double length = 0.0;         // path length inside the scintillator [mm]
    double Edep = 0.0;           // [MeV]
    double Evis = 0.0;           // [MeV]
    bool exits = false;          // left the scintillator before stopping
    std::int32_t n_hadronic = 0; // its own hadronic interactions (diagnostic)
    std::int32_t n_nuclear_recoil = 0; // its steps with DepositStep::nuclear_recoil
};

struct DepositRecord {
    double Edep = 0.0; // all particles [MeV]
    double Evis = 0.0; // all particles [MeV]
    ProtonTrack proton;
};

class DepositRecorder {
public:
    // Starts a new event, discarding the previous one.
    void beginEvent();
    void onStep(const DepositStep& step);

    [[nodiscard]] const DepositRecord& record() const {
        return record_;
    }

private:
    DepositRecord record_;
};

} // namespace riptide::sim
