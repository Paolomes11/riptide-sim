#include "riptide/sim/stepping_action.hpp"

#include <G4EmCalculator.hh>
#include <G4EmSaturation.hh>
#include <G4Exception.hh>
#include <G4HadronicProcess.hh>
#include <G4HadronicProcessType.hh>
#include <G4IonisParamMat.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Neutron.hh>
#include <G4Nucleus.hh>
#include <G4Proton.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4StepStatus.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserTrackInformation.hh>

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "riptide/common/units.hpp"
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

// Energy grid of the proton light table (D21): from well below the end of
// any proton range to five times the highest source energy of the project
// (20 MeV), 50 points per decade, so S changes by under 5% between points.
constexpr double kLightGridLow = 1.0 * units::keV;
constexpr double kLightGridHigh = 100.0 * units::MeV;
constexpr std::size_t kLightGridPointsPerDecade = 50;

// L(K) of protons in the scintillator from the mean stopping power of the
// proton's ionisation tables, the same tables that drive its transport.
// GetDEDX reads the built tables, so this runs after the physics is built;
// unlike ComputeNuclearDEDX it changes no model setting. Nuclear stopping,
// which only acts below 1 MeV in this physics list, is left out (D21).
ProtonLightTable buildProtonLightTable(const G4Material& material) {
    const double decades = std::log10(kLightGridHigh / kLightGridLow);
    const auto points =
        static_cast<std::size_t>(std::lround(decades * kLightGridPointsPerDecade)) + 1;
    std::vector<double> K(points);
    std::vector<double> S(points);
    G4EmCalculator calculator;
    for (std::size_t i = 0; i < points; ++i) {
        K[i] = kLightGridLow *
               std::pow(10.0, decades * static_cast<double>(i) / static_cast<double>(points - 1));
        // Project units equal Geant4's internal units (checked in sim/tests).
        S[i] = calculator.GetDEDX(K[i], G4Proton::Definition(), &material);
    }
    return {std::move(K), std::move(S), material.GetIonisation()->GetBirksConstant()};
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

double SteppingAction::visibleEnergy(const G4Step* step, const G4EmSaturation& saturation) {
    const G4Track* track = step->GetTrack();
    if (track->GetDefinition() != G4Proton::Definition()) {
        return saturation.VisibleEnergyDepositionAtAStep(step);
    }
    try {
        if (!proton_light_) {
            // Only steps starting in the scintillator get here: this is its material.
            proton_light_ = buildProtonLightTable(*step->GetPreStepPoint()->GetMaterial());
        }
        return proton_light_->visibleEnergy(step->GetPreStepPoint()->GetKineticEnergy(),
                                            step->GetTotalEnergyDeposit());
    } catch (const ProtonLightError& e) {
        // An empty dE/dx table or a proton above the grid: the light of the
        // run would be wrong, so stop it.
        G4Exception("riptide::sim::SteppingAction", "riptide_sim_002", FatalException,
                    (std::string("proton light: ") + e.what()).c_str());
        return 0.0;
    }
}

void SteppingAction::recordDeposit(const G4Step* step, bool hadronic) {
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
    d.evis = visibleEnergy(step, *saturation);
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
