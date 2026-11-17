#include "riptide/sim/stepping_action.hpp"

#include <G4HadronicProcess.hh>
#include <G4HadronicProcessType.hh>
#include <G4Nucleus.hh>
#include <G4Proton.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>

#include "riptide/sim/detector_construction.hpp"

namespace riptide::sim {

namespace {

// Geant4 numbers the primary particle of the event as track 1.
constexpr int kPrimaryTrackId = 1;

HadronicKind hadronicKind(int subtype) {
    switch (subtype) {
    case fHadronElastic:
        return HadronicKind::Elastic;
    case fHadronInelastic:
        return HadronicKind::Inelastic;
    case fCapture:
        return HadronicKind::Capture;
    case fFission:
        return HadronicKind::Fission;
    default:
        return HadronicKind::Other;
    }
}

Vec3 toVec3(const G4ThreeVector& v) {
    // Project units equal Geant4's internal units (checked in sim/tests).
    return {v.x(), v.y(), v.z()};
}

} // namespace

void SteppingAction::UserSteppingAction(const G4Step* step) {
    if (step->GetTrack()->GetTrackID() != kPrimaryTrackId) {
        return;
    }
    // The pre-step volume is where the step, and so the interaction, happened.
    const auto* volume = step->GetPreStepPoint()->GetPhysicalVolume();
    if (volume == nullptr || volume->GetName() != kScintillatorName) {
        return;
    }
    const auto* post = step->GetPostStepPoint();
    const auto* process = dynamic_cast<const G4HadronicProcess*>(post->GetProcessDefinedStep());
    if (process == nullptr) {
        return;
    }

    HadronicStep hadronic;
    hadronic.kind = hadronicKind(process->GetProcessSubType());
    // The process keeps the nucleus it sampled for this interaction until the next one.
    const auto* nucleus = process->GetTargetNucleus();
    hadronic.target_Z = nucleus != nullptr ? nucleus->GetZ_asInt() : 0;
    hadronic.position = toVec3(post->GetPosition());
    hadronic.time = post->GetGlobalTime();
    for (const auto* secondary : *step->GetSecondaryInCurrentStep()) {
        if (secondary->GetDefinition() == G4Proton::Definition()) {
            hadronic.proton = RecoilProton{secondary->GetKineticEnergy(),
                                           toVec3(secondary->GetMomentumDirection())};
            break;
        }
    }
    recorder_->onHadronicStep(hadronic);
}

} // namespace riptide::sim
