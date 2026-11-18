#include "riptide/sim/check_accumulator.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"

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

} // namespace

Nuclide makeNuclide(Target target, std::string label, double number_density,
                    const endf::EndfMaterial& material) {
    return Nuclide{target,
                   std::move(label),
                   number_density,
                   endf::crossSection(material, kTotalMt),
                   endf::crossSection(material, kElasticMt),
                   endf::crossSection(material, kCaptureMt)};
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

CheckSetup makeCheckSetup(const RunConfig& config, std::size_t bins_u, std::size_t bins_z) {
    if (bins_u == 0 || bins_z == 0) {
        throw std::invalid_argument("makeCheckSetup: bin counts must be positive");
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
    u_edges_ = uniformEdges(0.0, 1.0, setup_.bins_u);
    z_edges_ = uniformEdges(0.0, setup_.chord.length(), setup_.bins_z);
    const std::size_t n_models = models_.size();
    category_expected_.assign(n_models, {});
    u_observed_.assign(setup_.bins_u, 0.0);
    u_expected_.assign(n_models, std::vector<double>(setup_.bins_u, 0.0));
    z_observed_.assign(setup_.bins_z, 0.0);
    z_expected_.assign(n_models, std::vector<double>(setup_.bins_z, 0.0));
    coefficient_sums_.assign(n_models, {});
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
    constexpr double kHalfPi = 0.5 * units::pi;
    if (record.theta_p > kHalfPi) {
        ++counters_.theta_p_beyond_half_pi;
    }
    const double Ep_max = kinematics::maxProtonEnergy(record.En);
    const double u = record.Ep / Ep_max;
    const bool in_range = u <= 1.0;
    if (in_range) {
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
    // the Ep residual measures. protonEnergy <= Ep_max holds exactly, so the
    // clamp only removes rounding at theta_p = 0.
    const double Ep_kinematic = kinematics::protonEnergy(record.En, record.theta_p);
    mu_.push_back(std::clamp(1.0 - 2.0 * Ep_kinematic / Ep_max, -1.0, 1.0));
    Ep_residuals_.push_back(record.Ep - Ep_kinematic);
}

CheckResults CheckAccumulator::finish() const {
    CheckResults results;
    results.events = events_;
    results.np_events = mu_.size();
    results.counters = counters_;
    results.spectrum_u = compareBinned(u_edges_, u_observed_, models_, u_expected_);
    results.depth_z = compareBinned(z_edges_, z_observed_, models_, z_expected_);
    results.Ep_residuals = summarize(Ep_residuals_);

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
