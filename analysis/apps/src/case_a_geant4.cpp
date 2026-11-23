#include "riptide/apps/case_a_geant4.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "riptide/stats/split.hpp"

namespace riptide::apps {

namespace {

constexpr std::string_view kTruthOption = "--truth";
constexpr char kSeparator = ':';

std::vector<std::string> splitFields(std::string_view spec) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    for (;;) {
        const std::size_t end = spec.find(kSeparator, start);
        fields.emplace_back(spec.substr(start, end - start));
        if (end == std::string_view::npos) {
            return fields;
        }
        start = end + 1;
    }
}

double positiveEnergy(const std::string& text) {
    const auto value = parseCliNumber<double>(text, kTruthOption);
    if (!std::isfinite(value) || value <= 0.0) {
        throw CliError(std::string{kTruthOption} + ": energy '" + text +
                       "' is not positive and finite");
    }
    return value;
}

} // namespace

std::size_t pilotScaled(std::size_t count, bool pilot) noexcept {
    return pilot ? std::max(count / kPilotDivisor, kPilotMinEntries) : count;
}

stats::PseudoSample pseudoSampleFromEvents(const EventSample& events) {
    const std::size_t n = events.size();
    for (const auto* column : {&events.truth_En, &events.truth_Ep_phys, &events.truth_theta_p,
                               &events.reco_Ep, &events.reco_theta_p}) {
        if (column->size() != n) {
            throw std::invalid_argument(
                "pseudoSampleFromEvents: a column does not have one entry per event");
        }
    }
    stats::PseudoSample sample;
    sample.truth_En = events.truth_En;
    sample.truth_Ep = events.truth_Ep_phys;
    sample.truth_theta_p = events.truth_theta_p;
    sample.reco_Ep = events.reco_Ep;
    sample.reco_theta_p = events.reco_theta_p;
    return sample;
}

io::SourceSpectrum parseTruth(std::string_view spec) {
    const auto fields = splitFields(spec);
    io::SourceSpectrum truth;
    if (fields.size() == 2 && fields[0] == io::spectrumModeName(io::SpectrumMode::Mono)) {
        truth.mode = io::SpectrumMode::Mono;
        truth.En_MeV = positiveEnergy(fields[1]);
        return truth;
    }
    if (fields.size() == 3 && fields[0] == io::spectrumModeName(io::SpectrumMode::Flat)) {
        truth.mode = io::SpectrumMode::Flat;
        truth.En_min_MeV = positiveEnergy(fields[1]);
        truth.En_max_MeV = positiveEnergy(fields[2]);
        if (truth.En_min_MeV >= truth.En_max_MeV) {
            throw CliError(std::string{kTruthOption} + ": the flat range '" + std::string{spec} +
                           "' is empty");
        }
        return truth;
    }
    throw CliError(std::string{kTruthOption} + ": expected mono:<En> or flat:<min>:<max>, got '" +
                   std::string{spec} + "'");
}

std::size_t rowExperiments(std::size_t planned, std::size_t available_events,
                           std::size_t events_per_experiment, bool pilot) {
    const std::size_t fit = stats::maxDisjointExperiments(available_events, events_per_experiment);
    const std::size_t experiments = std::min(pilotScaled(planned, pilot), fit);
    return experiments < kMinExperiments ? 0 : experiments;
}

void requireMonoTruth(const stats::PseudoSample& sample, double En) {
    for (const double value : sample.truth_En) {
        if (!(std::abs(value - En) <= kMonoTruthRelTolerance * std::abs(En))) {
            throw std::invalid_argument("requireMonoTruth: an event has truth_En " +
                                        std::to_string(value) + " MeV, not the configured " +
                                        std::to_string(En) + " MeV");
        }
    }
}

} // namespace riptide::apps
