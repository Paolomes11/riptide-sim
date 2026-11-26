#pragma once

// Pseudo-experiments of the Case B model (D43), to validate the port against
// the toy model before any Geant4 data: one shared neutron direction per
// experiment, E_n drawn per event, single n-p scattering isotropic in the CM
// frame with exact kinematics, Gaussian noise on Ep and a Gaussian offset of
// the track in its tangent plane, as in the detector response. Event k draws
// from its own engine perEventEngine(seed, first_event_id + k), so an event
// does not depend on how many others are generated with it.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/common/random.hpp"
#include "riptide/stats/acceptance.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/grid.hpp"

namespace riptide::stats {

struct CaseBPseudoEvent {
    double truth_En = 0.0;
    double truth_Ep = 0.0;
    double truth_theta_p = 0.0; // angle between track and omega_n, rad
    geometry::Direction truth_track;
    double reco_Ep = 0.0;
    geometry::Direction reco_track;
};

// Uniform on the unit sphere: cos theta uniform in [-1, 1], then phi uniform
// in [0, 2 pi), in this order.
[[nodiscard]] geometry::Direction isotropicDirection(random::Engine& engine);

// One event at a given E_n with neutrons along the unit vector omega_n.
// Draws, in this order: cos theta_CM uniform in [-1, 1], hence Ep and theta_p
// from the exact kinematics; the azimuth of the track about omega_n, uniform;
// the noise on Ep; the track offset (smearDirection). Any E_n spectrum (Case
// C) reuses the kinematics and the detector response of Case B through it.
[[nodiscard]] CaseBPseudoEvent drawPseudoEventAtEnergy(double En,
                                                       const geometry::Direction& omega_n,
                                                       const CaseBResolution& resolution,
                                                       random::Engine& engine);

// The same with the selection acceptance of D40 at that E_n: the kinematics
// is drawn as above and kept with probability r(E_n, theta_p) / max r,
// repeating on the same engine until one is kept; the noise is drawn for the
// kept event only. The E_n of the event is the one given, so any E_n
// spectrum passes through unchanged.
[[nodiscard]] CaseBPseudoEvent drawPseudoEventAtEnergy(double En,
                                                       const geometry::Direction& omega_n,
                                                       const CaseBResolution& resolution,
                                                       const AcceptanceRatio& acceptance,
                                                       random::Engine& engine);

// One event with neutrons along the unit vector omega_n: E_n uniform on the
// domain, then drawPseudoEventAtEnergy on the same engine.
[[nodiscard]] CaseBPseudoEvent drawCaseBPseudoEvent(const EnergyDomain& domain,
                                                    const geometry::Direction& omega_n,
                                                    const CaseBResolution& resolution,
                                                    random::Engine& engine);

// The same with the selection acceptance of D40, as the toy twin of a
// Geant4 sample analysed with --acceptance (D45 note): E_n is drawn once,
// then drawPseudoEventAtEnergy with the acceptance on the same engine.
[[nodiscard]] CaseBPseudoEvent drawCaseBPseudoEvent(const EnergyDomain& domain,
                                                    const geometry::Direction& omega_n,
                                                    const CaseBResolution& resolution,
                                                    const AcceptanceRatio& acceptance,
                                                    random::Engine& engine);

struct CaseBPseudoSample {
    std::vector<double> truth_En;
    std::vector<double> truth_Ep;
    std::vector<double> truth_theta_p;
    std::vector<geometry::Direction> truth_track;
    std::vector<double> reco_Ep;
    std::vector<geometry::Direction> reco_track;

    [[nodiscard]] std::size_t size() const noexcept {
        return truth_En.size();
    }
};

// n events sharing omega_n. Throws std::invalid_argument for non-positive
// resolutions or a domain without 0 < En_min < En_max.
[[nodiscard]] CaseBPseudoSample caseBPseudoSample(const EnergyDomain& domain,
                                                  const geometry::Direction& omega_n, std::size_t n,
                                                  const CaseBResolution& resolution,
                                                  std::uint64_t seed,
                                                  std::uint64_t first_event_id = 0);

// The same with the selection acceptance (drawCaseBPseudoEvent above).
[[nodiscard]] CaseBPseudoSample
caseBPseudoSample(const EnergyDomain& domain, const geometry::Direction& omega_n, std::size_t n,
                  const CaseBResolution& resolution, std::uint64_t seed,
                  const AcceptanceRatio& acceptance, std::uint64_t first_event_id = 0);

struct CaseBPseudoExperiment {
    geometry::Direction truth_omega_n;
    CaseBPseudoSample sample;
};

// The truth omega_n isotropic, drawn by the engine of the id reserved to the
// truth (the last id of the block of experimentFirstEventId), then n events
// on the ids of the block. n must not exceed kMaxExperimentEvents.
[[nodiscard]] CaseBPseudoExperiment caseBPseudoExperiment(const EnergyDomain& domain, std::size_t n,
                                                          const CaseBResolution& resolution,
                                                          std::uint64_t seed, std::uint64_t section,
                                                          std::uint64_t index);

// The same at a known omega_n, as the toy twin of a fixed-direction Geant4
// run (D45).
[[nodiscard]] CaseBPseudoExperiment
caseBPseudoExperimentAtDirection(const geometry::Direction& omega_n, const EnergyDomain& domain,
                                 std::size_t n, const CaseBResolution& resolution,
                                 std::uint64_t seed, std::uint64_t section, std::uint64_t index);

} // namespace riptide::stats
