#pragma once

// Pseudo-experiments of the Case C model (D47): per experiment a truth
// (mu_E, sigma_E, omega_n), then per event E_n from the Gaussian spectrum and
// the kinematics and detector response of Case B at that E_n. The truth uses
// the engine of the id reserved to it, event k the k-th id of the block of
// experimentFirstEventId, as in Case B, so results do not depend on threads.

#include <cstddef>
#include <cstdint>

#include "riptide/common/direction.hpp"
#include "riptide/common/random.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"

namespace riptide::stats {

// Ranges of the truth draws: mu_E uniform, sigma_E log-uniform (MeV).
struct CaseCTruthRange {
    double mu_min = 0.0;
    double mu_max = 0.0;
    double sigma_min = 0.0;
    double sigma_max = 0.0;
};

// The checklist of D47, the same as the toy model.
inline constexpr CaseCTruthRange kCaseCChecklistTruth{13.0, 17.0, 0.8, 2.4};
// The edge row of D47: much of the spectrum beyond En_max = 20 MeV.
inline constexpr CaseCTruthRange kCaseCEdgeTruth{17.5, 19.0, 0.8, 2.4};

// How the per-event E_n follows the Gaussian spectrum.
enum class CaseCSpectrum {
    truncated,   // N(mu_E, sigma_E) truncated to the domain: the density of the model
    untruncated, // N(mu_E, sigma_E) as is: events outside the domain stay
};

struct CaseCTruth {
    double mu = 0.0;    // MeV
    double sigma = 0.0; // MeV
    geometry::Direction omega_n;
};

// Draws, in this order: mu_E, log sigma_E, omega_n (isotropicDirection).
// Throws std::invalid_argument unless 0 < sigma_min <= sigma_max and
// mu_min <= mu_max, all finite.
[[nodiscard]] CaseCTruth drawCaseCTruth(const CaseCTruthRange& range, random::Engine& engine);

// Smallest mass of N(mu, sigma) on the domain that truncated sampling
// accepts: below it the rejection loop would take more than 1/mass tries
// per event on average, a sign of a truth far outside the domain.
inline constexpr double kMinTruncatedMass = 1e-3;

// Mass of N(mu, sigma) on [En_min, En_max].
[[nodiscard]] double normalMassOnDomain(const EnergyDomain& domain, double mu, double sigma);

// E_n from N(mu, sigma), redrawn until it falls in [En_min, En_max]: exact
// for the truncated density. Throws std::invalid_argument for sigma <= 0 or
// a mass on the domain below kMinTruncatedMass.
[[nodiscard]] double drawTruncatedNormalEnergy(const EnergyDomain& domain, double mu, double sigma,
                                               random::Engine& engine);

// One event at the truth: E_n from the spectrum, then drawPseudoEventAtEnergy
// on the same engine. An untruncated E_n outside the domain is kept, but it
// must be positive for the kinematics: std::invalid_argument otherwise.
[[nodiscard]] CaseBPseudoEvent drawCaseCPseudoEvent(const EnergyDomain& domain,
                                                    const CaseCTruth& truth, CaseCSpectrum spectrum,
                                                    const CaseBResolution& resolution,
                                                    random::Engine& engine);

struct CaseCPseudoExperiment {
    CaseCTruth truth;
    CaseBPseudoSample sample;
};

// The truth from the engine of the id reserved to it, then n events on the
// ids of the block. n must not exceed kMaxExperimentEvents; throws
// std::invalid_argument for non-positive resolutions or a bad domain.
[[nodiscard]] CaseCPseudoExperiment
caseCPseudoExperiment(const EnergyDomain& domain, std::size_t n, const CaseBResolution& resolution,
                      const CaseCTruthRange& range, CaseCSpectrum spectrum, std::uint64_t seed,
                      std::uint64_t section, std::uint64_t index);

// The truth of experiment (section, index) as caseCPseudoExperiment draws
// it, with omega_n then replaced by the given direction: a Geant4 run has
// one fixed direction (D48), and its experiments and their toy twins share
// mu_E and sigma_E through this one function.
[[nodiscard]] CaseCTruth caseCTruthAtDirection(const CaseCTruthRange& range,
                                               const geometry::Direction& omega_n,
                                               std::uint64_t seed, std::uint64_t section,
                                               std::uint64_t index);

// The toy twin of a Geant4 experiment (D48): n events on the ids of the
// block of (section, index), E_n from the spectrum of the truth truncated to
// the domain, then drawPseudoEventAtEnergy. With an acceptance (non-null)
// the kinematics follows r of D40, as the selected Geant4 events do; E_n is
// drawn before it, so the detected spectrum stays the truncated Gaussian.
// Throws as caseCPseudoExperiment.
[[nodiscard]] CaseCPseudoExperiment
caseCPseudoExperimentAtTruth(const EnergyDomain& domain, std::size_t n,
                             const CaseBResolution& resolution, const CaseCTruth& truth,
                             std::uint64_t seed, std::uint64_t section, std::uint64_t index,
                             const AcceptanceRatio* acceptance = nullptr);

} // namespace riptide::stats
