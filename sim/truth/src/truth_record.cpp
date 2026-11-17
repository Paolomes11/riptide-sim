#include "riptide/sim/truth_record.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace riptide::sim {

using schema::FirstInteraction;

FirstInteraction classifyFirstInteraction(HadronicKind kind, int target_Z) {
    if (kind == HadronicKind::Elastic && target_Z == kHydrogenZ) {
        return FirstInteraction::NpElastic;
    }
    if (kind == HadronicKind::Elastic && target_Z == kCarbonZ) {
        return FirstInteraction::NCElastic;
    }
    if (kind == HadronicKind::Inelastic && target_Z == kCarbonZ) {
        return FirstInteraction::NCInelastic;
    }
    return FirstInteraction::Other;
}

double thetaP(const Vec3& omega_n, const Vec3& p_dir) {
    const double dot = omega_n.x * p_dir.x + omega_n.y * p_dir.y + omega_n.z * p_dir.z;
    const double norms =
        std::sqrt((omega_n.x * omega_n.x + omega_n.y * omega_n.y + omega_n.z * omega_n.z) *
                  (p_dir.x * p_dir.x + p_dir.y * p_dir.y + p_dir.z * p_dir.z));
    return std::acos(std::clamp(dot / norms, -1.0, 1.0));
}

void TruthRecorder::beginEvent(const EventMeta& meta, double En, const Vec3& omega_n) {
    meta_ = meta;
    record_ = TruthRecord{};
    record_.En = En;
    record_.omega_n = omega_n;
}

void TruthRecorder::onHadronicStep(const HadronicStep& step) {
    if (record_.first_interaction == FirstInteraction::Unknown) {
        record_.first_interaction = classifyFirstInteraction(step.kind, step.target_Z);
        record_.r0 = step.position;
        record_.t0 = step.time;
    }
    if (step.kind != HadronicKind::Elastic || step.target_Z != kHydrogenZ) {
        return;
    }
    ++record_.n_np_scatters;
    if (!step.proton) {
        ++record_.n_np_without_proton;
        return;
    }
    // Only the first n-p scatter defines the recoil proton; a missing proton on
    // that scatter leaves valid_np false rather than promoting a later one.
    if (record_.n_np_scatters == 1) {
        record_.Ep = step.proton->Ep;
        record_.p_dir = step.proton->dir;
        record_.theta_p = thetaP(record_.omega_n, step.proton->dir);
        record_.valid_np = true;
    }
}

void TruthTally::add(const TruthRecord& record) {
    ++events_;
    ++counts_.at(static_cast<std::size_t>(record.first_interaction));
    np_without_proton_ += record.n_np_without_proton;
}

} // namespace riptide::sim
