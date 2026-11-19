#include "riptide/sim/check_accumulator.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/hp_recoil.hpp"

namespace riptide::sim::checks {

namespace {

constexpr int kElasticMt = 2;
constexpr int kTotalMt = 1;
constexpr int kCaptureMt = 102;

// Two bins with kMinExpectedCount each are the least a chi^2 with one
// normalization constraint can use (ndf >= 1).
constexpr std::size_t kMinChi2Bins = 2;

std::vector<double> uniformEdges(double low, double high, std::size_t bins) {
    std::vector<double> edges(bins + 1);
    for (std::size_t i = 0; i <= bins; ++i) {
        edges.at(i) = low + (high - low) * static_cast<double>(i) / static_cast<double>(bins);
    }
    edges.back() = high; // exact upper edge, whatever the rounding of the step
    return edges;
}

// Bin of x among `bins` equal bins of [0, high], the upper edge in the last bin.
std::size_t binIndex(double x, double high, std::size_t bins) {
    const auto i = static_cast<std::size_t>(std::floor(x / high * static_cast<double>(bins)));
    return std::min(i, bins - 1);
}

void addTo(std::vector<double>& sums, const std::vector<double>& values) {
    if (sums.size() < values.size()) {
        sums.resize(values.size(), 0.0);
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
        sums.at(i) += values.at(i);
    }
}

bool sameDirection(const Vec3& a, const Vec3& b) {
    return std::abs(a.x - b.x) <= kDirectionTolerance &&
           std::abs(a.y - b.y) <= kDirectionTolerance && std::abs(a.z - b.z) <= kDirectionTolerance;
}

HistogramResult compareBinned(const std::vector<double>& edges, const std::vector<double>& observed,
                              const std::vector<ReferenceModel>& models,
                              const std::vector<std::vector<double>>& expected) {
    HistogramResult result{edges, observed, {}};
    for (std::size_t m = 0; m < models.size(); ++m) {
        HistogramCheck check{models.at(m).library, expected.at(m), std::nullopt, std::nullopt, {}};
        double total = 0.0;
        for (const double e : check.expected) {
            total += e;
        }
        if (total < kMinExpectedCount * static_cast<double>(kMinChi2Bins)) {
            check.skipped = "fewer expected events than two chi^2 bins need";
        } else {
            check.merged = mergeSparseBins(observed, check.expected, edges, kMinExpectedCount);
            if (check.merged->observed.size() < kMinChi2Bins) {
                check.skipped = "fewer than two bins after merging the sparse ones";
            } else {
                // The expected counts sum to the observed total by construction.
                check.comparison =
                    compareHistogram(check.merged->observed, check.merged->expected, 1);
            }
        }
        result.per_library.push_back(std::move(check));
    }
    return result;
}

FractionCheck
compareCategories(const std::string& library,
                  const std::array<std::uint64_t, kFirstInteractionCategories>& counts,
                  const std::array<double, kFirstInteractionCategories>& expected) {
    FractionCheck check;
    check.library = library;
    check.expected = expected;
    const auto absorbing = static_cast<std::size_t>(
        std::max_element(expected.begin(), expected.end()) - expected.begin());
    check.absorbing_category = static_cast<std::int32_t>(absorbing);

    std::vector<std::uint64_t> merged_counts;
    std::vector<double> merged_expected;
    std::uint64_t folded_counts = 0;
    double folded_expected = 0.0;
    for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
        if (c != absorbing && expected.at(c) < kMinExpectedCount) {
            folded_counts += counts.at(c);
            folded_expected += expected.at(c);
        }
    }
    double expected_total = 0.0;
    for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
        if (c == absorbing || expected.at(c) >= kMinExpectedCount) {
            const bool absorbs = c == absorbing;
            check.categories.push_back(static_cast<std::int32_t>(c));
            merged_counts.push_back(counts.at(c) + (absorbs ? folded_counts : 0));
            merged_expected.push_back(expected.at(c) + (absorbs ? folded_expected : 0.0));
            expected_total += merged_expected.back();
        }
    }
    if (check.categories.size() < kMinChi2Bins) {
        check.skipped = "fewer than two categories with enough expected events";
        return check;
    }
    // The expected counts add up to the number of events up to rounding;
    // dividing by their own sum makes the probabilities sum to 1 exactly.
    std::vector<double> probabilities;
    probabilities.reserve(merged_expected.size());
    for (const double e : merged_expected) {
        probabilities.push_back(e / expected_total);
    }
    check.comparison = compareFractions(merged_counts, probabilities);
    return check;
}

ResidualSummary summarize(const std::vector<double>& residuals) {
    ResidualSummary summary;
    summary.samples = residuals.size();
    if (residuals.empty()) {
        return summary;
    }
    double sum = 0.0;
    double sum_squares = 0.0;
    for (const double r : residuals) {
        sum += r;
        sum_squares += r * r;
        summary.max_abs = std::max(summary.max_abs, std::abs(r));
    }
    const auto n = static_cast<double>(residuals.size());
    summary.mean = sum / n;
    summary.rms = std::sqrt(sum_squares / n);
    return summary;
}

OffsetSummary summarizeOffset(const std::vector<double>& offsets,
                              const std::vector<double>& relative) {
    OffsetSummary summary;
    summary.absolute = summarize(offsets);
    if (relative.empty()) {
        return summary;
    }
    double sum = 0.0;
    for (const double r : relative) {
        sum += r;
    }
    summary.relative_mean = sum / static_cast<double>(relative.size());
    const auto [low, high] = std::minmax_element(relative.begin(), relative.end());
    summary.relative_min = *low;
    summary.relative_max = *high;
    return summary;
}

} // namespace

double nonHpElasticProbability(double En) noexcept {
    return std::clamp(
        (En - kChipsElasticMinEnergy) / (kHpElasticMaxEnergy - kChipsElasticMinEnergy), 0.0, 1.0);
}

Nuclide makeNuclide(Target target, std::string label, double number_density,
                    const endf::EndfMaterial& material) {
    return Nuclide{target,
                   std::move(label),
                   number_density,
                   endf::crossSection(material, kTotalMt),
                   endf::crossSection(material, kElasticMt),
                   endf::crossSection(material, kCaptureMt)};
}

namespace {

// The HP tables keep a discontinuity as a doubled energy, as ENDF does; the
// Tab1 grid accepts the same (validate rejects a third equal energy).
endf::Tab1 linLinTable(const HpCrossSection& table, const std::string& where) {
    endf::Tab1 tab;
    tab.boundaries = {static_cast<int>(table.En.size())};
    tab.laws = {endf::InterpolationLaw::LinLin};
    tab.x = table.En;
    tab.y = table.sigma;
    endf::validate(tab, where);
    return tab;
}

// Limit from below at a point of the grid: at a doubled energy the first of
// the two values, elsewhere the value itself (endf::evaluate gives the one
// from above).
double valueBelow(const endf::Tab1& tab, double x) {
    const auto it = std::lower_bound(tab.x.begin(), tab.x.end(), x);
    if (it != tab.x.end() && !(x < *it)) {
        return tab.y[static_cast<std::size_t>(std::distance(tab.x.begin(), it))];
    }
    return endf::evaluate(tab, x);
}

} // namespace

Nuclide makeHpNuclide(Target target, std::string label, double number_density,
                      const HpCrossSection& elastic, const HpCrossSection& inelastic,
                      const HpCrossSection& capture) {
    endf::Tab1 el = linLinTable(elastic, label + " elastic");
    const endf::Tab1 inel = linLinTable(inelastic, label + " inelastic");
    endf::Tab1 cap = linLinTable(capture, label + " capture");
    const std::initializer_list<const endf::Tab1*> parts{&el, &inel, &cap};
    const double low = std::max({el.x.front(), inel.x.front(), cap.x.front()});
    const double high = std::min({el.x.back(), inel.x.back(), cap.x.back()});
    if (!(low < high)) {
        throw HpDataError(label + ": the three cross sections share no energy interval");
    }
    std::vector<double> grid{low, high};
    for (const endf::Tab1* tab : parts) {
        for (const double x : tab->x) {
            if (x > low && x < high) {
                grid.push_back(x);
            }
        }
    }
    std::sort(grid.begin(), grid.end());
    endf::Tab1 total;
    for (const double x : grid) {
        // Each energy once, though several tables may list it.
        if (!total.x.empty() && !(total.x.back() < x)) {
            continue;
        }
        double below = 0.0;
        double above = 0.0;
        for (const endf::Tab1* tab : parts) {
            below += valueBelow(*tab, x);
            above += endf::evaluate(*tab, x);
        }
        // A jump of any part is a jump of the total: doubled energy, below
        // then above. Not at the ends, where only one side is tabulated.
        const bool jump = (below < above || above < below) && x > low && x < high;
        if (jump) {
            total.x.push_back(x);
            total.y.push_back(below);
        }
        total.x.push_back(x);
        total.y.push_back(above);
    }
    total.boundaries = {static_cast<int>(total.x.size())};
    total.laws = {endf::InterpolationLaw::LinLin};
    endf::validate(total, label + " total");
    return Nuclide{target,           std::move(label), number_density,
                   std::move(total), std::move(el),    std::move(cap)};
}

std::vector<double> ReferenceModel::coefficientsAt(double En) const {
    return endf::coefficients(np_angular, En);
}

std::vector<MediumComponent> ReferenceModel::mediumAt(double En) const {
    std::vector<MediumComponent> medium;
    medium.reserve(nuclides.size());
    for (const auto& nuclide : nuclides) {
        medium.push_back(MediumComponent{
            nuclide.target, nuclide.number_density, endf::evaluate(nuclide.elastic, En),
            endf::evaluate(nuclide.capture, En), endf::evaluate(nuclide.total, En)});
    }
    return medium;
}

CheckSetup makeCheckSetup(const RunConfig& config, double hp_target_mass, std::size_t bins_u,
                          std::size_t bins_z) {
    if (bins_u == 0 || bins_z == 0) {
        throw std::invalid_argument("makeCheckSetup: bin counts must be positive");
    }
    if (!std::isfinite(hp_target_mass) || !(hp_target_mass > 0.0)) {
        throw std::invalid_argument("makeCheckSetup: the HP target mass must be positive");
    }
    const Vec3 half{0.5 * config.geometry.scint_size.x, 0.5 * config.geometry.scint_size.y,
                    0.5 * config.geometry.scint_size.z};
    const auto chord = rayBoxChord(config.source.position, config.source.direction, half);
    if (!chord) {
        throw std::invalid_argument("makeCheckSetup: the beam misses the scintillator");
    }
    return CheckSetup{config.source.position,
                      config.source.direction,
                      half,
                      config.source.En_min,
                      config.source.En_max,
                      *chord,
                      hp_target_mass,
                      bins_u,
                      bins_z};
}

CheckAccumulator::CheckAccumulator(CheckSetup setup, std::vector<ReferenceModel> models)
    : setup_(setup), models_(std::move(models)) {
    if (models_.empty()) {
        throw std::invalid_argument("CheckAccumulator: no reference model");
    }
    std::set<std::string> names;
    for (const auto& model : models_) {
        if (!names.insert(model.library).second) {
            throw std::invalid_argument("CheckAccumulator: library '" + model.library +
                                        "' given twice");
        }
    }
    if (setup_.bins_u == 0 || setup_.bins_z == 0) {
        throw std::invalid_argument("CheckAccumulator: bin counts must be positive");
    }
    if (!std::isfinite(setup_.hp_target_mass) || !(setup_.hp_target_mass > 0.0)) {
        throw std::invalid_argument("CheckAccumulator: the HP target mass must be positive");
    }
    u_edges_ = uniformEdges(0.0, 1.0, setup_.bins_u);
    z_edges_ = uniformEdges(0.0, setup_.chord.length(), setup_.bins_z);
    const std::size_t n_models = models_.size();
    category_expected_.assign(n_models, {});
    u_observed_.assign(setup_.bins_u, 0.0);
    u_expected_.assign(n_models, std::vector<double>(setup_.bins_u, 0.0));
    z_observed_.assign(setup_.bins_z, 0.0);
    z_expected_.assign(n_models, std::vector<double>(setup_.bins_z, 0.0));
    coefficient_sums_.assign(n_models, {});

    const double mix_low = std::max(setup_.En_min, kChipsElasticMinEnergy);
    const double mix_high = std::min(setup_.En_max, kHpElasticMaxEnergy);
    if (mix_high > mix_low) {
        mix_edges_ = uniformEdges(mix_low, mix_high, kModelMixBins);
        mix_.classified_per_bin.assign(kModelMixBins, 0.0);
        mix_.observed_per_bin.assign(kModelMixBins, 0.0);
        mix_.expected_per_bin.assign(kModelMixBins, 0.0);
        mix_.variance_per_bin.assign(kModelMixBins, 0.0);
    }
}

void CheckAccumulator::add(const TruthRecord& record) {
    ++events_;
    const double tolerance = kEnergyRelTolerance * setup_.En_max;
    if (!(record.En >= setup_.En_min - tolerance && record.En <= setup_.En_max + tolerance)) {
        // The evaluated data would be read at an energy the run never asked
        // for; the event is counted and left out of every comparison.
        ++counters_.En_outside_source;
        return;
    }
    if (!sameDirection(record.omega_n, setup_.direction)) {
        ++counters_.omega_not_beam;
    }
    if (record.valid_np && record.n_np_scatters < 1) {
        ++counters_.valid_np_without_scatter;
    }

    std::vector<double> sigma_totals;
    sigma_totals.reserve(models_.size());
    for (std::size_t m = 0; m < models_.size(); ++m) {
        const auto probabilities =
            firstInteraction(models_.at(m).mediumAt(record.En), setup_.chord.length());
        for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
            category_expected_.at(m).at(c) += probabilities.unconditional.at(c);
        }
        sigma_totals.push_back(probabilities.sigma_total);
    }
    const auto category = static_cast<std::size_t>(record.first_interaction);
    ++category_counts_.at(category);

    if (record.first_interaction != schema::FirstInteraction::Unknown) {
        addDepth(record, sigma_totals);
    }
    if (record.first_interaction == schema::FirstInteraction::NpElastic) {
        if (record.valid_np) {
            addNpEvent(record);
        } else {
            ++counters_.np_first_without_proton;
        }
    }
}

void CheckAccumulator::addDepth(const TruthRecord& record,
                                const std::vector<double>& sigma_totals) {
    if (!insideBox(record.r0, setup_.half_size, kVertexTolerance)) {
        ++counters_.vertex_outside_box;
        return;
    }
    const double L = setup_.chord.length();
    // A vertex within the tolerance outside a face belongs to the edge bin.
    const double depth = std::clamp(
        depthAlongBeam(record.r0, setup_.position, setup_.direction, setup_.chord), 0.0, L);
    z_observed_.at(binIndex(depth, L, setup_.bins_z)) += 1.0;
    for (std::size_t m = 0; m < models_.size(); ++m) {
        addTo(z_expected_.at(m),
              truncatedExponentialBinProbabilities(sigma_totals.at(m), L, z_edges_));
    }
}

void CheckAccumulator::addNpEvent(const TruthRecord& record) {
    if (!std::isfinite(record.Ep) || !(record.Ep > 0.0)) {
        ++counters_.invalid_proton_energy;
        return;
    }
    // With the target at rest the proton cannot go backwards.
    constexpr double kHalfPi = 0.5 * units::pi;
    if (record.theta_p > kHalfPi) {
        ++counters_.theta_p_beyond_half_pi;
    }
    ++np_events_;
    const double tolerance = kModelRelTolerance * record.En;
    const double Ep_hp = hpProtonEnergy(record.En, record.theta_p, setup_.hp_target_mass);
    const double Ep_physical = kinematics::protonEnergy(record.En, record.theta_p);
    const bool hp = std::abs(record.Ep - Ep_hp) <= tolerance;
    const bool physical = std::abs(record.Ep - Ep_physical) <= tolerance;
    if (hp && physical) {
        // Too close to Ep = 0 to tell the models apart. Such events are too
        // few to matter and carry no information on the model, so they stay
        // in the HP sample and out of the model mix.
        ++counters_.np_model_ambiguous;
        addHpEvent(record);
        return;
    }
    if (!hp && !physical) {
        ++counters_.np_kinematics_unmatched;
        return;
    }
    addModelMix(record.En, physical);
    if (physical) {
        if (!(nonHpElasticProbability(record.En) > 0.0)) {
            ++counters_.np_non_hp_outside_window;
        }
        Ep_non_hp_residuals_.push_back(record.Ep - Ep_physical);
        return;
    }
    addHpEvent(record);
}

void CheckAccumulator::addModelMix(double En, bool non_hp) {
    const double p = nonHpElasticProbability(En);
    ++mix_.classified;
    if (non_hp) {
        ++mix_.observed_non_hp;
    }
    mix_.expected_non_hp += p;
    mix_.variance += p * (1.0 - p);
    if (mix_edges_.empty() || En < mix_edges_.front()) {
        return;
    }
    const double low = mix_edges_.front();
    const std::size_t bin = binIndex(En - low, mix_edges_.back() - low, kModelMixBins);
    mix_.classified_per_bin.at(bin) += 1.0;
    mix_.observed_per_bin.at(bin) += non_hp ? 1.0 : 0.0;
    mix_.expected_per_bin.at(bin) += p;
    mix_.variance_per_bin.at(bin) += p * (1.0 - p);
}

// Every check below compares the HP final state with the evaluated data, so
// only HP events enter it. Leaving out the non-HP ones keeps the expectations
// exact: the model is drawn from E_n before the angle, and the expectations
// are summed at each kept event's own E_n.
void CheckAccumulator::addHpEvent(const TruthRecord& record) {
    const double M_T = setup_.hp_target_mass;
    // u is the CM variable of the evaluated angular data, recovered exactly
    // from Ep through the HP kinematics; it is in (0, 1] up to rounding.
    const double u_raw = 0.5 * (1.0 - hpCosCm(record.Ep, record.En, M_T));
    const bool in_range = u_raw <= 1.0 + kUnitTolerance;
    if (in_range) {
        const double u = std::min(u_raw, 1.0);
        u_observed_.at(binIndex(u, 1.0, setup_.bins_u)) += 1.0;
    } else {
        ++counters_.u_overflow;
    }
    for (std::size_t m = 0; m < models_.size(); ++m) {
        const auto a = models_.at(m).coefficientsAt(record.En);
        addTo(coefficient_sums_.at(m), a);
        if (in_range) {
            // Only histogrammed events add to the expectation, so observed and
            // expected totals agree and the chi^2 has one constraint.
            addTo(u_expected_.at(m), protonEnergyBinProbabilities(a, u_edges_));
        }
    }
    // mu from the angle, not from Ep: the difference between the two is what
    // the Ep residual measures. The clamp only removes rounding at theta_p = 0.
    mu_.push_back(std::clamp(hpCosCmFromAngle(record.En, record.theta_p, M_T), -1.0, 1.0));
    Ep_residuals_.push_back(record.Ep - hpProtonEnergy(record.En, record.theta_p, M_T));
    // The departure from the physical n-p kinematics, kept for the record:
    // it is the HP target-mass convention, not an error of the simulation.
    const double Ep_physical = kinematics::protonEnergy(record.En, record.theta_p);
    Ep_physical_offsets_.push_back(record.Ep - Ep_physical);
    if (Ep_physical > 0.0) {
        Ep_physical_relative_.push_back(record.Ep / Ep_physical - 1.0);
    }
}

ModelMixResult CheckAccumulator::finishModelMix() const {
    ModelMixResult mix = mix_;
    mix.edges = mix_edges_;
    mix.Ep_residuals = summarize(Ep_non_hp_residuals_);
    // The normal approximation of a sum of Bernoulli trials needs a variance
    // of a few counts, as the chi^2 bins need kMinExpectedCount.
    if (mix.variance < kMinExpectedCount) {
        mix.skipped = "binomial variance of the non-HP count below the minimum expected count";
    } else {
        mix.pull = (static_cast<double>(mix.observed_non_hp) - mix.expected_non_hp) /
                   std::sqrt(mix.variance);
    }
    if (mix_edges_.empty()) {
        mix.chi2_skipped = "the source range does not overlap the HP/CHIPS window";
        return mix;
    }
    double chi2 = 0.0;
    for (std::size_t i = 0; i < kModelMixBins; ++i) {
        const double variance = mix.variance_per_bin.at(i);
        if (variance < kMinExpectedCount) {
            mix.chi2_skipped =
                "a window bin has a binomial variance below the minimum expected count";
            return mix;
        }
        const double d = mix.observed_per_bin.at(i) - mix.expected_per_bin.at(i);
        chi2 += d * d / variance;
    }
    // The bins are independent sums of trials and nothing is fitted: ndf = bins.
    mix.ndf = static_cast<int>(kModelMixBins);
    mix.chi2 = chi2;
    mix.p_value = chi2SurvivalProbability(chi2, mix.ndf);
    return mix;
}

CheckResults CheckAccumulator::finish() const {
    CheckResults results;
    results.events = events_;
    results.np_events = np_events_;
    results.np_hp_events = mu_.size();
    results.counters = counters_;
    results.model_mix = finishModelMix();
    results.spectrum_u = compareBinned(u_edges_, u_observed_, models_, u_expected_);
    results.depth_z = compareBinned(z_edges_, z_observed_, models_, z_expected_);
    results.Ep_residuals = summarize(Ep_residuals_);
    results.Ep_physical_offset = summarizeOffset(Ep_physical_offsets_, Ep_physical_relative_);

    results.fractions.observed = category_counts_;
    std::uint64_t accepted = 0;
    for (const auto count : category_counts_) {
        accepted += count;
    }
    for (std::size_t m = 0; m < models_.size(); ++m) {
        if (accepted == 0) {
            FractionCheck check;
            check.library = models_.at(m).library;
            check.skipped = "no events";
            results.fractions.per_library.push_back(std::move(check));
        } else {
            results.fractions.per_library.push_back(compareCategories(
                models_.at(m).library, category_counts_, category_expected_.at(m)));
        }
    }

    std::size_t max_degree = 0;
    for (const auto& sums : coefficient_sums_) {
        max_degree = std::max(max_degree, sums.size());
    }
    if (mu_.size() < 2) {
        results.legendre.skipped = "fewer than two n-p events";
        return results;
    }
    if (max_degree == 0) {
        results.legendre.skipped = "the reference angular distributions have no coefficients";
        return results;
    }
    auto moments = legendreMoments(mu_, static_cast<int>(max_degree));
    for (const double error : moments.errors) {
        if (!(error > 0.0)) {
            results.legendre.skipped = "a Legendre moment has zero spread in the sample";
            return results;
        }
    }
    const auto n = static_cast<double>(mu_.size());
    for (std::size_t m = 0; m < models_.size(); ++m) {
        LegendreCheck check;
        check.library = models_.at(m).library;
        check.expected.assign(max_degree, 0.0); // a_l beyond the file's order are zero
        for (std::size_t l = 0; l < coefficient_sums_.at(m).size(); ++l) {
            check.expected.at(l) = coefficient_sums_.at(m).at(l) / n;
        }
        for (std::size_t l = 0; l < max_degree; ++l) {
            check.pulls.push_back((moments.values.at(l) - check.expected.at(l)) /
                                  moments.errors.at(l));
        }
        results.legendre.per_library.push_back(std::move(check));
    }
    results.legendre.moments = std::move(moments);
    return results;
}

} // namespace riptide::sim::checks
