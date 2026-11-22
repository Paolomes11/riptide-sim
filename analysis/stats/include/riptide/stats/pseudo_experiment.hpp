#pragma once

// Pseudo-experiments of the Case A model, to validate the statistics engine
// against itself before any Geant4 data: single n-p scattering, isotropic in
// the CM frame, exact kinematics, and 1D Gaussian noise on Ep and theta_p as
// in the detector response. Event k draws from its own engine
// perEventEngine(seed, first_event_id + k), so an event does not depend on
// how many others are generated with it.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "riptide/common/random.hpp"
#include "riptide/stats/case_a.hpp"
#include "riptide/stats/grid.hpp"

namespace riptide::stats {

struct PseudoEvent {
    double truth_En = 0.0;
    double truth_Ep = 0.0;
    double truth_theta_p = 0.0;
    double reco_Ep = 0.0;
    double reco_theta_p = 0.0;
};

// One event at a known E_n: cos theta_CM uniform in [-1, 1], hence Ep
// uniform in [0, Ep_max(E_n)), then the recoil angle from the exact inverse
// kinematics and the measurement noise, in this order.
[[nodiscard]] PseudoEvent drawPseudoEvent(double En, const CaseAResolution& resolution,
                                          random::Engine& engine);

struct PseudoSample {
    std::vector<double> truth_En;
    std::vector<double> truth_Ep;
    std::vector<double> truth_theta_p;
    std::vector<double> reco_Ep;
    std::vector<double> reco_theta_p;

    [[nodiscard]] std::size_t size() const noexcept {
        return truth_En.size();
    }
};

// n events with E_n drawn per event from the flat prior on the domain, as in
// a calibration with the truth taken from the prior.
[[nodiscard]] PseudoSample pseudoSampleFromPrior(const EnergyDomain& domain, std::size_t n,
                                                 const CaseAResolution& resolution,
                                                 std::uint64_t seed,
                                                 std::uint64_t first_event_id = 0);

// n events sharing one E_n, as in a combined measurement.
[[nodiscard]] PseudoSample pseudoSampleAtEnergy(double En, std::size_t n,
                                                const CaseAResolution& resolution,
                                                std::uint64_t seed,
                                                std::uint64_t first_event_id = 0);

} // namespace riptide::stats
