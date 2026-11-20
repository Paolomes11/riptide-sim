#include "riptide/sim/pstar_checks.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

#include "riptide/sim/binning.hpp"

namespace riptide::sim::checks {

void ResidualStats::add(double observed, double expected) {
    ++samples_;
    const double n = static_cast<double>(samples_);
    observed_mean_ += (observed - observed_mean_) / n;
    expected_mean_ += (expected - expected_mean_) / n;
    const double residual = observed - expected;
    const double delta = residual - residual_mean_;
    residual_mean_ += delta / n;
    residual_m2_ += delta * (residual - residual_mean_);
}

std::optional<double> ResidualStats::residualSd() const {
    if (samples_ < 2) {
        return std::nullopt;
    }
    return std::sqrt(residual_m2_ / static_cast<double>(samples_ - 1));
}

std::optional<double> ResidualStats::residualSem() const {
    const std::optional<double> sd = residualSd();
    if (!sd) {
        return std::nullopt;
    }
    return *sd / std::sqrt(static_cast<double>(samples_));
}

std::optional<double> ResidualStats::relativeResidual() const {
    if (samples_ == 0 || !(expected_mean_ > 0.0)) {
        return std::nullopt;
    }
    return residual_mean_ / expected_mean_;
}

std::optional<double> ResidualStats::relativeSem() const {
    const std::optional<double> sem = residualSem();
    if (!sem || !(expected_mean_ > 0.0)) {
        return std::nullopt;
    }
    return *sem / expected_mean_;
}

PstarComparison::PstarComparison(const DepositSetup& setup, PstarTable table)
    : En_max_(setup.En_max), bins_(setup.bins), light_(std::move(table), setup.birks_kB) {
    if (!(En_max_ > 0.0)) {
        throw std::invalid_argument("PstarComparison: non-positive En_max");
    }
    if (bins_ == 0) {
        throw std::invalid_argument("PstarComparison: zero bins");
    }
    if (!(En_max_ <= light_.table().Ep.back())) {
        throw PstarError("PstarComparison: En_max above the PSTAR table");
    }
    results_.edges = uniformEdges(0.0, En_max_, bins_);
    results_.range_bins.assign(bins_, {});
    results_.light_bins.assign(bins_, {});
}

void PstarComparison::add(const TruthRecord& truth, const DepositRecord& deposit) {
    const ProtonTrack& p = deposit.proton;
    // A nuclear recoil shortens the track and takes its energy off the light,
    // so such a proton is not the pure electronic slowing-down PSTAR describes.
    if (!truth.valid_np || p.exits || p.n_hadronic != 0 || p.n_nuclear_recoil != 0) {
        return;
    }
    const double Ep = truth.Ep;
    // The CSDA range is not tabulated below the first energy.
    if (!(Ep >= light_.table().Ep.front())) {
        ++results_.below_table;
        return;
    }
    const double range = csdaRange(light_.table(), Ep);
    const double light = light_.light(Ep);
    ++results_.samples;
    // The HP recoil can exceed En_max by its small offset: last bin, as in
    // the deposit checks.
    const std::size_t bin = binIndex(Ep, En_max_, bins_);
    results_.range_bins.at(bin).add(p.length, range);
    results_.light_bins.at(bin).add(p.Evis, light);
    results_.range.add(p.length, range);
    results_.light.add(p.Evis, light);
}

} // namespace riptide::sim::checks
