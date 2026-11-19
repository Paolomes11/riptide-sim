#include "riptide/sim/stepping_action.hpp"

#include <G4EmSaturation.hh>
#include <G4Exception.hh>
#include <G4HadronicProcess.hh>
#include <G4HadronicProcessType.hh>
#include <G4Neutron.hh>
#include <G4Nucleus.hh>
#include <G4Proton.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4StepStatus.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserTrackInformation.hh>

#include "riptide/sim/detector_construction.hpp"

namespace riptide::sim {

namespace {

// Geant4 numbers the primary particle of the event as track 1.
constexpr int kPrimaryTrackId = 1;

// Tag carried by the recoil proton of the first n-p scatter. The track owns
// it and deletes it with itself; secondaries of the proton do not inherit it.
class RecoilProtonInfo final : public G4VUserTrackInformation {};

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

TrackRole trackRole(const G4Track& track) {
    return dynamic_cast<const RecoilProtonInfo*>(track.GetUserInformation()) != nullptr
               ? TrackRole::RecoilProton
               : TrackRole::Other;
}

} // namespace

void SteppingAction::UserSteppingAction(const G4Step* step) {
    // The pre-step volume is where the step, and so any interaction, happened.
    if (step->GetPreStepPoint()->GetPhysicalVolume() != detector_->scintillator()) {
        return;
    }
    const auto* process =
        dynamic_cast<const G4HadronicProcess*>(step->GetPostStepPoint()->GetProcessDefinedStep());
    recordDeposit(step, process != nullptr);
    if (process == nullptr) {
        return;
    }
    const G4Track* track = step->GetTrack();
    if (track->GetTrackID() == kPrimaryTrackId) {
        onPrimaryHadronic(step, *process);
    } else if (track->GetDefinition() == G4Neutron::Definition()) {
        truth_->onSecondaryNeutronHadronic();
    }
}

void SteppingAction::recordDeposit(const G4Step* step, bool hadronic) const {
    const G4EmSaturation* saturation = detector_->emSaturation();
    if (saturation == nullptr) {
        G4Exception("riptide::sim::SteppingAction", "riptide_sim_001", FatalException,
                    "G4EmSaturation not created: DetectorConstruction::Construct has not run");
        return;
    }
    const auto* pre = step->GetPreStepPoint();
    const auto* post = step->GetPostStepPoint();
    DepositStep d;
    d.role = trackRole(*step->GetTrack());
    d.edep = step->GetTotalEnergyDeposit();
    d.evis = saturation->VisibleEnergyDepositionAtAStep(step);
    d.pre = toVec3(pre->GetPosition());
    d.post = toVec3(post->GetPosition());
    d.length = step->GetStepLength();
    // The scintillator has no daughter volumes, so a step that ends on a
    // geometric boundary ends on its surface and the track goes on outside.
    d.leaves_volume = post->GetStepStatus() == fGeomBoundary;
    d.hadronic = hadronic;
    deposits_->onStep(d);
}

void SteppingAction::onPrimaryHadronic(const G4Step* step, const G4HadronicProcess& process) const {
    const auto* post = step->GetPostStepPoint();
    HadronicStep hadronic;
    hadronic.kind = hadronicKind(process.GetProcessSubType());
    // The process keeps the nucleus it sampled for this interaction until the next one.
    const auto* nucleus = process.GetTargetNucleus();
    hadronic.target_Z = nucleus != nullptr ? nucleus->GetZ_asInt() : 0;
    hadronic.position = toVec3(post->GetPosition());
    hadronic.time = post->GetGlobalTime();
    const G4Track* proton = nullptr;
    for (const auto* secondary : *step->GetSecondaryInCurrentStep()) {
        if (secondary->GetDefinition() == G4Proton::Definition()) {
            proton = secondary;
            hadronic.proton = RecoilProton{secondary->GetKineticEnergy(),
                                           toVec3(secondary->GetMomentumDirection())};
            break;
        }
    }
    const bool had_np = truth_->record().valid_np;
    truth_->onHadronicStep(hadronic);
    // The recorder decides which n-p is "the" one (the first with a proton);
    // the proton it took its Ep from is the one to follow.
    if (!had_np && truth_->record().valid_np && proton != nullptr) {
        proton->SetUserInformation(new RecoilProtonInfo);
    }
}

} // namespace riptide::sim
