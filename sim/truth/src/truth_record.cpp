#include "riptide/sim/truth_record.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace riptide::sim {

using schema::EventClass;
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

InteractionCounts interactionCounts(const TruthRecord& record) {
    return {record.n_np_scatters, record.n_nC_elastic, record.n_nC_inelastic, record.n_other};
}

EventClass classifyEvent(const InteractionCounts& counts) {
    if (counts.np < 0 || counts.nC_elastic < 0 || counts.nC_inelastic < 0 || counts.other < 0) {
        throw std::invalid_argument("classifyEvent: negative interaction count");
    }
    const bool has_np = counts.np > 0;
    const bool has_nC = counts.nC_elastic > 0 || counts.nC_inelastic > 0;
    if (counts.other > 0) {
        return EventClass::Other;
    }
    if (has_np && has_nC) {
        return EventClass::Mixed;
    }
    if (has_np) {
        return counts.np == 1 ? EventClass::SingleNp : EventClass::MultipleNp;
    }
    if (counts.nC_inelastic > 0) {
        return EventClass::NCInelasticNoNp;
    }
    if (counts.nC_elastic > 0) {
        return EventClass::NCElasticOnly;
    }
    return EventClass::NoInteraction;
}

double thetaP(const Vec3& omega_n, const Vec3& p_dir) {
    const double dot = omega_n.x * p_dir.x + omega_n.y * p_dir.y + omega_n.z * p_dir.z;
    const double norms =
        std::sqrt((omega_n.x * omega_n.x + omega_n.y * omega_n.y + omega_n.z * omega_n.z) *
                  (p_dir.x * p_dir.x + p_dir.y * p_dir.y + p_dir.z * p_dir.z));
    return std::acos(std::clamp(dot / norms, -1.0, 1.0));
}

void TruthRecorder::beginEvent(const EventMeta& meta, double En, const Vec3& omega_n,
                               const Vec3& n_start) {
    meta_ = meta;
    record_ = TruthRecord{};
    record_.En = En;
    record_.omega_n = omega_n;
    record_.n_start = n_start;
}

void TruthRecorder::onHadronicStep(const HadronicStep& step) {
    if (record_.first_interaction == FirstInteraction::Unknown) {
        record_.first_interaction = classifyFirstInteraction(step.kind, step.target_Z);
        record_.r0 = step.position;
        record_.t0 = step.time;
    }
    switch (classifyFirstInteraction(step.kind, step.target_Z)) {
    case FirstInteraction::NpElastic:
        break;
    case FirstInteraction::NCElastic:
        ++record_.n_nC_elastic;
        return;
    case FirstInteraction::NCInelastic:
        ++record_.n_nC_inelastic;
        return;
    case FirstInteraction::Unknown:
    case FirstInteraction::Other:
        ++record_.n_other;
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
        record_.np_vertex = step.position;
        record_.valid_np = true;
    }
}

void TruthRecorder::onSecondaryNeutronHadronic() {
    ++record_.n_secondary_n_hadronic;
}

void TruthTally::add(const TruthRecord& record) {
    ++events_;
    ++counts_.at(static_cast<std::size_t>(record.first_interaction));
    np_without_proton_ += record.n_np_without_proton;
}

} // namespace riptide::sim
