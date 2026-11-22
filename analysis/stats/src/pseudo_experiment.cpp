#include "riptide/stats/pseudo_experiment.hpp"

#include <cmath>
#include <stdexcept>

#include "riptide/common/np_kinematics.hpp"

namespace riptide::stats {

namespace {

constexpr double kCosSpan = 2.0;

void requireValid(const CaseAResolution& resolution) {
    if (!(resolution.sigma_Ep > 0.0 && resolution.sigma_theta > 0.0)) {
        throw std::invalid_argument("pseudo-experiments need positive resolutions");
    }
}

void append(PseudoSample& sample, const PseudoEvent& event) {
    sample.truth_En.push_back(event.truth_En);
    sample.truth_Ep.push_back(event.truth_Ep);
    sample.truth_theta_p.push_back(event.truth_theta_p);
    sample.reco_Ep.push_back(event.reco_Ep);
    sample.reco_theta_p.push_back(event.reco_theta_p);
}

void reserve(PseudoSample& sample, std::size_t n) {
    sample.truth_En.reserve(n);
    sample.truth_Ep.reserve(n);
    sample.truth_theta_p.reserve(n);
    sample.reco_Ep.reserve(n);
    sample.reco_theta_p.reserve(n);
}

} // namespace

PseudoEvent drawPseudoEvent(double En, const CaseAResolution& resolution, random::Engine& engine) {
    if (!(En > 0.0) || !std::isfinite(En)) {
        throw std::invalid_argument("pseudo-events need a positive finite E_n");
    }
    PseudoEvent event;
    event.truth_En = En;
    const double cos_cm = 1.0 - (kCosSpan * random::uniform53(engine));
    event.truth_Ep = kinematics::protonEnergyFromCosCm(En, cos_cm);
    event.truth_theta_p = kinematics::protonAngleFromEnergy(event.truth_Ep, En);
    event.reco_Ep = event.truth_Ep + (resolution.sigma_Ep * random::standardNormal(engine));
    event.reco_theta_p =
        event.truth_theta_p + (resolution.sigma_theta * random::standardNormal(engine));
    return event;
}

PseudoSample pseudoSampleFromPrior(const EnergyDomain& domain, std::size_t n,
                                   const CaseAResolution& resolution, std::uint64_t seed,
                                   std::uint64_t first_event_id) {
    requireValid(resolution);
    if (!(domain.En_min > 0.0 && domain.En_max > domain.En_min)) {
        throw std::invalid_argument("pseudo-experiments need 0 < En_min < En_max");
    }
    PseudoSample sample;
    reserve(sample, n);
    for (std::size_t k = 0; k < n; ++k) {
        auto engine = random::perEventEngine(seed, first_event_id + k);
        const double En =
            domain.En_min + ((domain.En_max - domain.En_min) * random::uniform53(engine));
        append(sample, drawPseudoEvent(En, resolution, engine));
    }
    return sample;
}

PseudoSample pseudoSampleAtEnergy(double En, std::size_t n, const CaseAResolution& resolution,
                                  std::uint64_t seed, std::uint64_t first_event_id) {
    requireValid(resolution);
    PseudoSample sample;
    reserve(sample, n);
    for (std::size_t k = 0; k < n; ++k) {
        auto engine = random::perEventEngine(seed, first_event_id + k);
        append(sample, drawPseudoEvent(En, resolution, engine));
    }
    return sample;
}

} // namespace riptide::stats
