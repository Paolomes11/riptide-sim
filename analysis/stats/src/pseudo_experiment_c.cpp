#include "riptide/stats/pseudo_experiment_c.hpp"

#include <cmath>
#include <stdexcept>

#include "riptide/stats/calibration.hpp"

namespace riptide::stats {

namespace {

void requireValid(const EnergyDomain& domain, const CaseBResolution& resolution) {
    if (!(resolution.sigma_Ep > 0.0 && resolution.sigma_theta > 0.0)) {
        throw std::invalid_argument("Case C pseudo-experiments need positive resolutions");
    }
    if (!(domain.En_min > 0.0 && domain.En_max > domain.En_min) || !std::isfinite(domain.En_max)) {
        throw std::invalid_argument("Case C pseudo-experiments need 0 < En_min < En_max");
    }
}

void requireSpectrum(double mu, double sigma) {
    if (!std::isfinite(mu) || !std::isfinite(sigma) || !(sigma > 0.0)) {
        throw std::invalid_argument("the Case C spectrum needs a finite mu and sigma > 0");
    }
}

double uniformBetween(double lo, double hi, random::Engine& engine) {
    return lo + ((hi - lo) * random::uniform53(engine));
}

double drawNormalEnergy(double mu, double sigma, random::Engine& engine) {
    return mu + (sigma * random::standardNormal(engine));
}

void requireBlockSize(std::size_t n) {
    if (n > kMaxExperimentEvents) {
        throw std::invalid_argument("too many events for one pseudo-experiment block");
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

} // namespace

CaseCTruth drawCaseCTruth(const CaseCTruthRange& range, random::Engine& engine) {
    const bool finite = std::isfinite(range.mu_min) && std::isfinite(range.mu_max) &&
                        std::isfinite(range.sigma_min) && std::isfinite(range.sigma_max);
    if (!finite || !(range.mu_min <= range.mu_max) ||
        !(range.sigma_min > 0.0 && range.sigma_min <= range.sigma_max)) {
        throw std::invalid_argument("Case C truth range needs mu_min <= mu_max, "
                                    "0 < sigma_min <= sigma_max");
    }
    CaseCTruth truth;
    truth.mu = uniformBetween(range.mu_min, range.mu_max, engine);
    truth.sigma =
        std::exp(uniformBetween(std::log(range.sigma_min), std::log(range.sigma_max), engine));
    truth.omega_n = isotropicDirection(engine);
    return truth;
}

double normalMassOnDomain(const EnergyDomain& domain, double mu, double sigma) {
    requireSpectrum(mu, sigma);
    // The mass is Phi(z_hi) - Phi(z_lo), z = (E - mu) / sigma. A domain all
    // above mu is written as Q(z_lo) - Q(z_hi), one all below as
    // Phi(z_hi) - Phi(z_lo) with Phi(z) = erfc(-z / sqrt 2) / 2: erfc then
    // has positive arguments and keeps the relative precision of a far tail,
    // where Phi or 1 - Phi would round to 0 or 1.
    const double scale = sigma * std::sqrt(2.0);
    const double t_lo = (domain.En_min - mu) / scale;
    const double t_hi = (domain.En_max - mu) / scale;
    if (t_lo > 0.0) {
        return 0.5 * (std::erfc(t_lo) - std::erfc(t_hi));
    }
    return 0.5 * (std::erfc(-t_hi) - std::erfc(-t_lo));
}

double drawTruncatedNormalEnergy(const EnergyDomain& domain, double mu, double sigma,
                                 random::Engine& engine) {
    if (normalMassOnDomain(domain, mu, sigma) < kMinTruncatedMass) {
        throw std::invalid_argument("the Case C spectrum has too little mass on the domain");
    }
    while (true) {
        const double En = drawNormalEnergy(mu, sigma, engine);
        if (En >= domain.En_min && En <= domain.En_max) {
            return En;
        }
    }
}

CaseBPseudoEvent drawCaseCPseudoEvent(const EnergyDomain& domain, const CaseCTruth& truth,
                                      CaseCSpectrum spectrum, const CaseBResolution& resolution,
                                      random::Engine& engine) {
    double En = 0.0;
    switch (spectrum) {
    case CaseCSpectrum::truncated:
        En = drawTruncatedNormalEnergy(domain, truth.mu, truth.sigma, engine);
        break;
    case CaseCSpectrum::untruncated:
        requireSpectrum(truth.mu, truth.sigma);
        En = drawNormalEnergy(truth.mu, truth.sigma, engine);
        break;
    }
    if (!(En > 0.0)) {
        throw std::invalid_argument("an untruncated Case C event has E_n <= 0");
    }
    return drawPseudoEventAtEnergy(En, truth.omega_n, resolution, engine);
}

CaseCPseudoExperiment caseCPseudoExperiment(const EnergyDomain& domain, std::size_t n,
                                            const CaseBResolution& resolution,
                                            const CaseCTruthRange& range, CaseCSpectrum spectrum,
                                            std::uint64_t seed, std::uint64_t section,
                                            std::uint64_t index) {
    requireValid(domain, resolution);
    requireBlockSize(n);
    const std::uint64_t first = experimentFirstEventId(section, index);
    auto truth_engine = random::perEventEngine(seed, first + kMaxExperimentEvents);
    CaseCPseudoExperiment experiment;
    experiment.truth = drawCaseCTruth(range, truth_engine);
    reserve(experiment.sample, n);
    for (std::size_t k = 0; k < n; ++k) {
        auto engine = random::perEventEngine(seed, first + k);
        append(experiment.sample,
               drawCaseCPseudoEvent(domain, experiment.truth, spectrum, resolution, engine));
    }
    return experiment;
}

CaseCTruth caseCTruthAtDirection(const CaseCTruthRange& range, const geometry::Direction& omega_n,
                                 std::uint64_t seed, std::uint64_t section, std::uint64_t index) {
    const std::uint64_t first = experimentFirstEventId(section, index);
    auto engine = random::perEventEngine(seed, first + kMaxExperimentEvents);
    CaseCTruth truth = drawCaseCTruth(range, engine);
    truth.omega_n = omega_n;
    return truth;
}

CaseCPseudoExperiment caseCPseudoExperimentAtTruth(const EnergyDomain& domain, std::size_t n,
                                                   const CaseBResolution& resolution,
                                                   const CaseCTruth& truth, std::uint64_t seed,
                                                   std::uint64_t section, std::uint64_t index,
                                                   const AcceptanceRatio* acceptance) {
    requireValid(domain, resolution);
    requireBlockSize(n);
    const std::uint64_t first = experimentFirstEventId(section, index);
    CaseCPseudoExperiment experiment;
    experiment.truth = truth;
    reserve(experiment.sample, n);
    for (std::size_t k = 0; k < n; ++k) {
        auto engine = random::perEventEngine(seed, first + k);
        const double En = drawTruncatedNormalEnergy(domain, truth.mu, truth.sigma, engine);
        append(experiment.sample,
               acceptance == nullptr
                   ? drawPseudoEventAtEnergy(En, truth.omega_n, resolution, engine)
                   : drawPseudoEventAtEnergy(En, truth.omega_n, resolution, *acceptance, engine));
    }
    return experiment;
}

} // namespace riptide::stats
