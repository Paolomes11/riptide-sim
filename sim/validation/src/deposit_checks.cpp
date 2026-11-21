#include "riptide/sim/deposit_checks.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "riptide/sim/binning.hpp"

namespace riptide::sim::checks {

namespace {

using schema::FirstInteraction;

constexpr double kHalf = 0.5;

// Distance of p outside the box along the axis that sticks out most: negative
// inside (minus the distance to the nearest face), zero on a face.
double surfaceDistance(const Vec3& p, const Vec3& half) {
    return std::max({std::abs(p.x) - half.x, std::abs(p.y) - half.y, std::abs(p.z) - half.z});
}

// The counter of the kind of the first interaction is nonzero: the first
// interaction is also counted among all of them.
bool firstInteractionCounted(const TruthRecord& truth) {
    const InteractionCounts c = interactionCounts(truth);
    switch (truth.first_interaction) {
    case FirstInteraction::Unknown:
        return c.np == 0 && c.nC_elastic == 0 && c.nC_inelastic == 0 && c.other == 0;
    case FirstInteraction::NpElastic:
        return c.np > 0;
    case FirstInteraction::NCElastic:
        return c.nC_elastic > 0;
    case FirstInteraction::NCInelastic:
        return c.nC_inelastic > 0;
    case FirstInteraction::Other:
        return c.other > 0;
    }
    return false;
}

void addToBin(std::vector<std::uint64_t>& counts, std::uint64_t& overflow, double x, double high) {
    if (x <= 0.0) {
        return;
    }
    if (x > high) {
        ++overflow;
        return;
    }
    ++counts.at(binIndex(x, high, counts.size()));
}

} // namespace

DepositSetup makeDepositSetup(const RunConfig& config, std::size_t bins) {
    if (!config.scintillator) {
        throw std::invalid_argument(
            "makeDepositSetup: the configuration has no scintillator block (version 1)");
    }
    const Vec3& size = config.geometry.scint_size;
    return {{kHalf * size.x, kHalf * size.y, kHalf * size.z},
            config.source.En_max,
            config.scintillator->birks_kB,
            bins};
}

DepositAccumulator::DepositAccumulator(DepositSetup setup) : setup_(setup) {
    const Vec3& h = setup_.half_size;
    if (!(h.x > 0.0 && h.y > 0.0 && h.z > 0.0)) {
        throw std::invalid_argument("DepositAccumulator: non-positive scintillator edge");
    }
    if (!(setup_.En_max > 0.0)) {
        throw std::invalid_argument("DepositAccumulator: non-positive En_max");
    }
    if (!(setup_.birks_kB >= 0.0)) {
        throw std::invalid_argument("DepositAccumulator: negative Birks constant");
    }
    if (setup_.bins == 0) {
        throw std::invalid_argument("DepositAccumulator: zero bins");
    }
    results_.edges = uniformEdges(0.0, setup_.En_max, setup_.bins);
    results_.Edep_counts.assign(setup_.bins, 0);
    results_.Evis_counts.assign(setup_.bins, 0);
    results_.protons_per_Ep.assign(setup_.bins, 0);
    results_.exiting_per_Ep.assign(setup_.bins, 0);
}

bool DepositAccumulator::onSurface(const Vec3& p) const {
    return std::abs(surfaceDistance(p, setup_.half_size)) <= kSurfaceTolerance;
}

bool DepositAccumulator::insideBox(const Vec3& p) const {
    return surfaceDistance(p, setup_.half_size) <= kSurfaceTolerance;
}

void DepositAccumulator::add(const TruthRecord& truth, const DepositRecord& deposit,
                             std::int32_t stored_class) {
    DepositCounters& c = results_.counters;
    ++results_.events;

    const ProtonTrack& p = deposit.proton;
    if (std::min({deposit.Edep, deposit.Evis, p.Edep, p.Evis}) < -kDepositSumTolerance) {
        ++c.negative_deposit;
    }
    // Birks quenching only ever lowers the deposit of a step.
    if (deposit.Evis > deposit.Edep + kDepositSumTolerance) {
        ++c.Evis_above_Edep;
    }
    const double visible_difference = std::abs(deposit.Evis - deposit.Edep);
    results_.kB_zero_max_abs_difference =
        std::max(results_.kB_zero_max_abs_difference, visible_difference);
    // With kB = 0 G4EmSaturation returns each step's deposit unchanged, and
    // both totals add the same steps in the same order.
    if (setup_.birks_kB <= 0.0 && visible_difference > kDepositSumTolerance) {
        ++c.kB_zero_visible_mismatch;
    }

    const InteractionCounts counts = interactionCounts(truth);
    const schema::EventClass event_class = classifyEvent(counts);
    ++results_.classes.at(static_cast<std::size_t>(event_class));
    if (stored_class != static_cast<std::int32_t>(event_class)) {
        ++c.class_mismatch;
    }
    if (!firstInteractionCounted(truth)) {
        ++c.first_interaction_uncounted;
    }
    // A neutron that never interacted leaves no secondary to deposit anything.
    if (truth.first_interaction == FirstInteraction::Unknown && deposit.Edep > 0.0) {
        ++c.deposit_without_interaction;
    }
    if (truth.n_secondary_n_hadronic > 0) {
        ++results_.events_secondary_n_hadronic;
    }

    if (truth.valid_np) {
        addProton(truth, deposit);
    }
    addSpectra(deposit);
}

void DepositAccumulator::addProton(const TruthRecord& truth, const DepositRecord& deposit) {
    DepositCounters& c = results_.counters;
    const ProtonTrack& p = deposit.proton;
    ++results_.np_events;

    if (p.Evis > p.Edep + kDepositSumTolerance) {
        ++c.p_Evis_above_p_Edep;
    }
    if (p.Edep > deposit.Edep + kDepositSumTolerance) {
        ++c.p_Edep_above_det_Edep;
    }
    const double excess = p.Edep - truth.Ep;
    if (excess > kNuclearStoppingSlack) {
        ++c.p_Edep_above_Ep_slack;
    } else if (excess > kDepositSumTolerance) {
        ++c.p_Edep_above_Ep;
    }
    if (!(p.length > 0.0)) {
        ++c.proton_without_track;
    }
    if (!insideBox(p.end)) {
        ++c.end_outside_box;
    }
    if (p.exits && !onSurface(p.end)) {
        ++c.exit_not_on_surface;
    }
    if (!p.exits && onSurface(p.end)) {
        ++c.contained_end_on_surface;
    }

    if (p.n_hadronic > 0) {
        ++results_.protons_hadronic;
    }
    if (p.n_nuclear_recoil > 0) {
        ++results_.protons_nuclear_recoil;
    }
    const std::size_t bin = binIndex(std::max(truth.Ep, 0.0), setup_.En_max, setup_.bins);
    ++results_.protons_per_Ep.at(bin);
    if (p.exits) {
        ++results_.protons_exiting;
        ++results_.exiting_per_Ep.at(bin);
    }

    if (!p.exits && p.n_hadronic == 0 && p.n_nuclear_recoil == 0) {
        BalanceSummary& b = results_.contained_balance;
        b.min = b.samples == 0 ? excess : std::min(b.min, excess);
        b.max = b.samples == 0 ? excess : std::max(b.max, excess);
        ++b.samples;
        balance_sum_ += excess;
        if (excess > kDepositSumTolerance) {
            ++b.positive;
        } else if (excess < -kDepositSumTolerance) {
            ++b.negative;
        }
    }
}

void DepositAccumulator::addSpectra(const DepositRecord& deposit) {
    addToBin(results_.Edep_counts, results_.Edep_overflow, deposit.Edep, setup_.En_max);
    addToBin(results_.Evis_counts, results_.Evis_overflow, deposit.Evis, setup_.En_max);
}

DepositResults DepositAccumulator::finish() const {
    DepositResults results = results_;
    BalanceSummary& b = results.contained_balance;
    if (b.samples > 0) {
        b.mean = balance_sum_ / static_cast<double>(b.samples);
    }
    return results;
}

} // namespace riptide::sim::checks
