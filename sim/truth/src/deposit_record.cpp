#include "riptide/sim/deposit_record.hpp"

namespace riptide::sim {

void DepositRecorder::beginEvent() {
    record_ = DepositRecord{};
}

void DepositRecorder::onStep(const DepositStep& step) {
    record_.Edep += step.edep;
    record_.Evis += step.evis;
    if (step.role != TrackRole::RecoilProton) {
        return;
    }
    ProtonTrack& p = record_.proton;
    if (!p.seen) {
        p.seen = true;
        p.start = step.pre;
    }
    // Steps of one track arrive in order, so the last one seen gives the end
    // point and tells whether the proton left the volume. A proton cannot
    // come back: the block is convex and the world around it is vacuum, so a
    // proton that leaves moves on a straight line away from it.
    p.end = step.post;
    p.exits = step.leaves_volume;
    p.length += step.length;
    p.Edep += step.edep;
    p.Evis += step.evis;
    if (step.hadronic) {
        ++p.n_hadronic;
    }
    if (step.nuclear_recoil) {
        ++p.n_nuclear_recoil;
    }
}

} // namespace riptide::sim
