#include "riptide/stats/pseudo_experiment_b.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/calibration.hpp"

namespace riptide::stats {

namespace {

constexpr double kCosSpan = 2.0;

void requireValid(const EnergyDomain& domain, const CaseBResolution& resolution) {
    if (!(resolution.sigma_Ep > 0.0 && resolution.sigma_theta > 0.0)) {
        throw std::invalid_argument("Case B pseudo-experiments need positive resolutions");
    }
    if (!(domain.En_min > 0.0 && domain.En_max > domain.En_min) || !std::isfinite(domain.En_max)) {
        throw std::invalid_argument("Case B pseudo-experiments need 0 < En_min < En_max");
    }
}

void reserve(CaseBPseudoSample& sample, std::size_t n) {
    sample.truth_En.reserve(n);
    sample.truth_Ep.reserve(n);
    sample.truth_theta_p.reserve(n);
    sample.truth_track.reserve(n);
    sample.reco_Ep.reserve(n);
    sample.reco_track.reserve(n);
}

void append(CaseBPseudoSample& sample, const CaseBPseudoEvent& event) {
    sample.truth_En.push_back(event.truth_En);
    sample.truth_Ep.push_back(event.truth_Ep);
    sample.truth_theta_p.push_back(event.truth_theta_p);
    sample.truth_track.push_back(event.truth_track);
    sample.reco_Ep.push_back(event.reco_Ep);
    sample.reco_track.push_back(event.reco_track);
}

void requireBlockSize(std::size_t n) {
    if (n > kMaxExperimentEvents) {
        throw std::invalid_argument("too many events for one pseudo-experiment block");
    }
}

} // namespace

geometry::Direction isotropicDirection(random::Engine& engine) {
    const double cos_theta = 1.0 - (kCosSpan * random::uniform53(engine));
    const double phi = 2.0 * units::pi * random::uniform53(engine);
    // max() guards 1 - cos^2 against a rounding just below zero at the poles.
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - (cos_theta * cos_theta)));
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
}

CaseBPseudoEvent drawCaseBPseudoEvent(const EnergyDomain& domain,
                                      const geometry::Direction& omega_n,
                                      const CaseBResolution& resolution, random::Engine& engine) {
    CaseBPseudoEvent event;
    event.truth_En = domain.En_min + ((domain.En_max - domain.En_min) * random::uniform53(engine));
    const double cos_cm = 1.0 - (kCosSpan * random::uniform53(engine));
    event.truth_Ep = kinematics::protonEnergyFromCosCm(event.truth_En, cos_cm);
    event.truth_theta_p = kinematics::protonAngleFromEnergy(event.truth_Ep, event.truth_En);
    const double phi = 2.0 * units::pi * random::uniform53(engine);
    event.truth_track = geometry::directionAroundAxis(omega_n, event.truth_theta_p, phi);
    event.reco_Ep = event.truth_Ep + (resolution.sigma_Ep * random::standardNormal(engine));
    event.reco_track = geometry::smearDirection(event.truth_track, resolution.sigma_theta, engine);
    return event;
}

CaseBPseudoSample caseBPseudoSample(const EnergyDomain& domain, const geometry::Direction& omega_n,
                                    std::size_t n, const CaseBResolution& resolution,
                                    std::uint64_t seed, std::uint64_t first_event_id) {
    requireValid(domain, resolution);
    CaseBPseudoSample sample;
    reserve(sample, n);
    for (std::size_t k = 0; k < n; ++k) {
        auto engine = random::perEventEngine(seed, first_event_id + k);
        append(sample, drawCaseBPseudoEvent(domain, omega_n, resolution, engine));
    }
    return sample;
}

CaseBPseudoExperiment caseBPseudoExperiment(const EnergyDomain& domain, std::size_t n,
                                            const CaseBResolution& resolution, std::uint64_t seed,
                                            std::uint64_t section, std::uint64_t index) {
    requireBlockSize(n);
    const std::uint64_t first = experimentFirstEventId(section, index);
    auto engine = random::perEventEngine(seed, first + kMaxExperimentEvents);
    const geometry::Direction omega_n = isotropicDirection(engine);
    return {omega_n, caseBPseudoSample(domain, omega_n, n, resolution, seed, first)};
}

CaseBPseudoExperiment caseBPseudoExperimentAtDirection(const geometry::Direction& omega_n,
                                                       const EnergyDomain& domain, std::size_t n,
                                                       const CaseBResolution& resolution,
                                                       std::uint64_t seed, std::uint64_t section,
                                                       std::uint64_t index) {
    requireBlockSize(n);
    return {omega_n, caseBPseudoSample(domain, omega_n, n, resolution, seed,
                                       experimentFirstEventId(section, index))};
}

} // namespace riptide::stats
