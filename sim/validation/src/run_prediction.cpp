#include "riptide/sim/run_prediction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/check_report.hpp"
#include "riptide/sim/hp_recoil.hpp"

namespace riptide::sim::checks {

namespace {

// Newton iterations for the Gauss-Legendre nodes converge quadratically from
// the asymptotic guess; the loop stops at a step of a few ulp, and the cap
// only guards against a failure to converge.
constexpr int kNewtonMaxIterations = 100;
constexpr double kNewtonTolerance = 4.0 * std::numeric_limits<double>::epsilon();

// Bisection for the normal quantile: the bracket covers every p down to
// erfc's underflow, and 200 halvings reach the double resolution.
constexpr double kQuantileUpper = 40.0;
constexpr int kBisectionSteps = 200;

// Proton lab angles sampled for the range of the HP offset. The offset is
// smooth in theta_p, so the grid resolves its extremes far below the
// tolerance of any threshold built on them; the last point stays below pi/2,
// where Ep vanishes and the ratio is 0/0.
constexpr std::size_t kOffsetAngles = 4096;

constexpr double kHalfPi = 0.5 * units::pi;

std::size_t category(schema::FirstInteraction c) {
    return static_cast<std::size_t>(c);
}

std::vector<double> uniformEdges(double low, double high, std::size_t bins) {
    std::vector<double> edges(bins + 1);
    for (std::size_t i = 0; i <= bins; ++i) {
        edges.at(i) = low + (high - low) * static_cast<double>(i) / static_cast<double>(bins);
    }
    edges.back() = high;
    return edges;
}

void addScaled(std::vector<double>& sum, const std::vector<double>& values, double weight) {
    if (sum.size() < values.size()) {
        sum.resize(values.size(), 0.0);
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
        sum.at(i) += weight * values.at(i);
    }
}

// Legendre polynomial and its derivative at x for the Newton step.
std::pair<double, double> legendreWithDerivative(int degree, double x) {
    const double p = kinematics::legendreP(degree, x);
    const double previous = kinematics::legendreP(degree - 1, x);
    return {p, static_cast<double>(degree) * (x * p - previous) / (x * x - 1.0)};
}

void offsetRange(const CheckSetup& setup, const std::vector<double>& energies,
                 RunPrediction& run) {
    double low = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();
    for (const double En : energies) {
        for (std::size_t k = 0; k < kOffsetAngles; ++k) {
            const double theta_p =
                kHalfPi * static_cast<double>(k) / static_cast<double>(kOffsetAngles);
            const double physical = kinematics::protonEnergy(En, theta_p);
            const double hp = hpProtonEnergy(En, theta_p, setup.hp_target_mass);
            const double relative = (hp - physical) / physical;
            low = std::min(low, relative);
            high = std::max(high, relative);
        }
    }
    run.offset_relative_min = low;
    run.offset_relative_max = high;
}

LibraryPrediction predictLibrary(const CheckSetup& setup, const ReferenceModel& model,
                                 const QuadratureRule& spectrum) {
    const double L = setup.chord.length();
    const std::vector<double> u_edges = [&] {
        return uniformEdges(0.0, 1.0, setup.bins_u);
    }();
    LibraryPrediction p;
    p.library = model.library;
    p.u_probabilities.assign(setup.bins_u, 0.0);
    std::vector<double> second_moments;
    for (std::size_t i = 0; i < spectrum.nodes.size(); ++i) {
        const double En = spectrum.nodes.at(i);
        const double w = spectrum.weights.at(i);
        const auto first = firstInteraction(model.mediumAt(En), L);
        p.sigma_total += w * first.sigma_total;
        for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
            p.probabilities.at(c) += w * first.unconditional.at(c);
        }
        const double np = first.unconditional.at(category(schema::FirstInteraction::NpElastic));
        const double non_hp = nonHpElasticProbability(En);
        const double hp_weight = w * np * (1.0 - non_hp);
        p.hp_np_probability += hp_weight;
        p.non_hp_np_probability += w * np * non_hp;
        p.non_hp_variance += w * np * non_hp * (1.0 - non_hp);
        if (!(hp_weight > 0.0)) {
            continue;
        }
        const auto a = model.coefficientsAt(En);
        addScaled(p.legendre, a, hp_weight);
        std::vector<double> squares(a.size());
        for (std::size_t l = 0; l < a.size(); ++l) {
            squares.at(l) = legendreSecondMoment(a, static_cast<int>(l + 1));
        }
        addScaled(second_moments, squares, hp_weight);
        addScaled(p.u_probabilities, protonEnergyBinProbabilities(a, u_edges), hp_weight);
    }
    if (!(p.hp_np_probability > 0.0)) {
        throw std::invalid_argument("predictRun: the source gives no HP n-p events");
    }
    for (double& value : p.legendre) {
        value /= p.hp_np_probability;
    }
    for (double& value : p.u_probabilities) {
        value /= p.hp_np_probability;
    }
    p.legendre_spread.resize(p.legendre.size());
    for (std::size_t l = 0; l < p.legendre.size(); ++l) {
        const double mean = p.legendre.at(l);
        const double variance = second_moments.at(l) / p.hp_np_probability - mean * mean;
        p.legendre_spread.at(l) = std::sqrt(std::max(variance, 0.0));
    }
    return p;
}

double maxAbs(const std::vector<double>& values) {
    double m = 0.0;
    for (const double v : values) {
        m = std::max(m, std::abs(v));
    }
    return m;
}

nlohmann::json categoriesJson(const std::array<double, kFirstInteractionCategories>& values) {
    nlohmann::json doc = nlohmann::json::object();
    for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
        doc[std::string(categoryName(static_cast<std::int32_t>(c)))] = values.at(c);
    }
    return doc;
}

} // namespace

QuadratureRule gaussLegendre(std::size_t points) {
    if (points == 0) {
        throw std::invalid_argument("gaussLegendre: at least one point");
    }
    const int n = static_cast<int>(points);
    QuadratureRule rule;
    rule.nodes.resize(points);
    rule.weights.resize(points);
    for (int i = 0; i < n; ++i) {
        // Asymptotic guess of the i-th root, descending from near +1.
        double x = std::cos(units::pi * (i + 0.75) / (n + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < kNewtonMaxIterations; ++iteration) {
            const auto [p, dp] = legendreWithDerivative(n, x);
            derivative = dp;
            const double step = p / dp;
            x -= step;
            if (std::abs(step) <= kNewtonTolerance) {
                break;
            }
        }
        derivative = legendreWithDerivative(n, x).second;
        const auto index = static_cast<std::size_t>(n - 1 - i); // ascending nodes
        rule.nodes.at(index) = x;
        rule.weights.at(index) = 2.0 / ((1.0 - x * x) * derivative * derivative);
    }
    return rule;
}

QuadratureRule sourceQuadrature(double En_min, double En_max, std::size_t points) {
    if (!(En_min > 0.0) || En_max < En_min || points == 0) {
        throw std::invalid_argument("sourceQuadrature: need 0 < En_min <= En_max and points > 0");
    }
    QuadratureRule rule;
    if (!(En_max > En_min)) {
        rule.nodes = {En_min};
        rule.weights = {1.0};
        return rule;
    }
    const double step = (En_max - En_min) / static_cast<double>(points);
    rule.nodes.resize(points);
    rule.weights.assign(points, 1.0 / static_cast<double>(points));
    for (std::size_t i = 0; i < points; ++i) {
        rule.nodes.at(i) = En_min + (static_cast<double>(i) + 0.5) * step;
    }
    return rule;
}

double legendreSecondMoment(const std::vector<double>& a, int degree) {
    if (degree < 1) {
        throw std::invalid_argument("legendreSecondMoment: degree must be at least 1");
    }
    // P_degree^2 f has degree 2 degree + a.size(); n points integrate up to 2n - 1.
    const auto points = static_cast<std::size_t>(degree) + a.size() + 1;
    const QuadratureRule rule = gaussLegendre(points);
    double integral = 0.0;
    for (std::size_t i = 0; i < points; ++i) {
        const double mu = rule.nodes.at(i);
        double f = 0.5;
        for (std::size_t l = 0; l < a.size(); ++l) {
            const auto order = static_cast<int>(l + 1);
            f += 0.5 * (2.0 * order + 1.0) * a.at(l) * kinematics::legendreP(order, mu);
        }
        const double p = kinematics::legendreP(degree, mu);
        integral += rule.weights.at(i) * p * p * f;
    }
    return integral;
}

double normalTwoSidedQuantile(double p) {
    if (!(p > 0.0 && p < 1.0)) {
        throw std::invalid_argument("normalTwoSidedQuantile: p must lie in (0, 1)");
    }
    // P(|Z| > z) = erfc(z / sqrt 2), decreasing in z.
    double low = 0.0;
    double high = kQuantileUpper;
    for (int step = 0; step < kBisectionSteps; ++step) {
        const double mid = 0.5 * (low + high);
        if (std::erfc(mid / std::sqrt(2.0)) > p) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return 0.5 * (low + high);
}

RunPrediction predictRun(const CheckSetup& setup, const std::vector<ReferenceModel>& models,
                         std::uint64_t events, std::size_t energy_nodes) {
    if (models.empty() || events == 0 || energy_nodes == 0) {
        throw std::invalid_argument("predictRun: need models, events and energy nodes");
    }
    const QuadratureRule spectrum = sourceQuadrature(setup.En_min, setup.En_max, energy_nodes);
    RunPrediction run;
    run.events = events;
    run.L = setup.chord.length();
    run.En_min = setup.En_min;
    run.En_max = setup.En_max;
    run.energy_nodes = spectrum.nodes.size();
    run.overlaps_model_window =
        setup.En_max > kChipsElasticMinEnergy && setup.En_min < kHpElasticMaxEnergy;
    std::vector<double> energies = spectrum.nodes;
    energies.push_back(setup.En_min);
    energies.push_back(setup.En_max);
    offsetRange(setup, energies, run);
    std::size_t degrees = 0;
    for (const auto& model : models) {
        run.libraries.push_back(predictLibrary(setup, model, spectrum));
        degrees = std::max(degrees, run.libraries.back().legendre.size());
    }
    // As in the checks, a_l beyond a file's order are zero.
    for (auto& library : run.libraries) {
        library.legendre.resize(degrees, 0.0);
        library.legendre_spread.resize(degrees, 0.0);
    }
    return run;
}

ExpectedErrors expectedErrors(const RunPrediction& run, std::size_t reference) {
    const LibraryPrediction& lib = run.libraries.at(reference);
    const auto n = static_cast<double>(run.events);
    ExpectedErrors errors;
    errors.hp_np_events = n * lib.hp_np_probability;
    for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
        const double p = lib.probabilities.at(c);
        errors.fraction_sigma.at(c) = std::sqrt(p * (1.0 - p) / n);
    }
    for (const double spread : lib.legendre_spread) {
        errors.legendre_sigma.push_back(spread / std::sqrt(errors.hp_np_events));
    }
    return errors;
}

PredictionShift predictionShift(const LibraryPrediction& reference, const LibraryPrediction& other,
                                const ExpectedErrors& errors, std::uint64_t events) {
    if (reference.legendre.size() != other.legendre.size() ||
        reference.legendre.size() != errors.legendre_sigma.size() ||
        reference.u_probabilities.size() != other.u_probabilities.size()) {
        throw std::invalid_argument("predictionShift: predictions of different shape");
    }
    PredictionShift shift;
    std::vector<double> scaled;
    for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
        const double sigma = errors.fraction_sigma.at(c);
        if (sigma > 0.0) {
            scaled.push_back((other.probabilities.at(c) - reference.probabilities.at(c)) / sigma);
        }
    }
    shift.fractions = maxAbs(scaled);
    scaled.clear();
    for (std::size_t l = 0; l < reference.legendre.size(); ++l) {
        scaled.push_back((other.legendre.at(l) - reference.legendre.at(l)) /
                         errors.legendre_sigma.at(l));
    }
    shift.legendre = maxAbs(scaled);
    const double hp_events = static_cast<double>(events) * reference.hp_np_probability;
    for (std::size_t b = 0; b < reference.u_probabilities.size(); ++b) {
        const double q = reference.u_probabilities.at(b);
        if (q > 0.0) {
            const double d = other.u_probabilities.at(b) - q;
            shift.u_noncentrality += hp_events * d * d / q;
        }
    }
    return shift;
}

TestPlan testPlan(const RunPrediction& run, double alpha_family) {
    if (!(alpha_family > 0.0 && alpha_family < 1.0)) {
        throw std::invalid_argument("testPlan: alpha must lie in (0, 1)");
    }
    TestPlan plan;
    constexpr int kHistogramAndFractionTests = 3; // u, depth, fractions
    plan.chi2_tests = kHistogramAndFractionTests + (run.overlaps_model_window ? 1 : 0);
    const auto degrees =
        run.libraries.empty() ? 0 : static_cast<int>(run.libraries.front().legendre.size());
    plan.pull_tests = degrees + (run.overlaps_model_window ? 1 : 0);
    plan.alpha_family = alpha_family;
    plan.alpha_per_test = alpha_family / (plan.chi2_tests + plan.pull_tests);
    plan.pull_threshold = normalTwoSidedQuantile(plan.alpha_per_test);
    return plan;
}

nlohmann::json predictionJson(const RunPrediction& run, const PredictionShift& quadrature,
                              const TestPlan& plan) {
    const auto n = static_cast<double>(run.events);
    nlohmann::json doc;
    doc["events"] = run.events;
    doc["L_mm"] = run.L;
    doc["En_min_MeV"] = run.En_min;
    doc["En_max_MeV"] = run.En_max;
    doc["energy_nodes"] = run.energy_nodes;
    doc["overlaps_model_window"] = run.overlaps_model_window;
    doc["Ep_physical_offset_relative"] = {{"min", run.offset_relative_min},
                                          {"max", run.offset_relative_max}};
    doc["quadrature_shift_in_sigma"] = {{"fractions", quadrature.fractions},
                                        {"legendre", quadrature.legendre},
                                        {"u_noncentrality", quadrature.u_noncentrality}};
    doc["test_plan"] = {{"chi2_tests", plan.chi2_tests},
                        {"pull_tests", plan.pull_tests},
                        {"alpha_family", plan.alpha_family},
                        {"alpha_per_test", plan.alpha_per_test},
                        {"pull_threshold", plan.pull_threshold}};
    doc["libraries"] = nlohmann::json::array();
    for (std::size_t m = 0; m < run.libraries.size(); ++m) {
        const LibraryPrediction& lib = run.libraries.at(m);
        const ExpectedErrors errors = expectedErrors(run, m);
        std::array<double, kFirstInteractionCategories> counts{};
        std::array<double, kFirstInteractionCategories> relative{};
        for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
            counts.at(c) = n * lib.probabilities.at(c);
            const double p = lib.probabilities.at(c);
            relative.at(c) = p > 0.0 ? errors.fraction_sigma.at(c) / p : 0.0;
        }
        nlohmann::json entry;
        entry["library"] = lib.library;
        entry["sigma_total_per_mm"] = lib.sigma_total;
        entry["probabilities"] = categoriesJson(lib.probabilities);
        entry["expected_counts"] = categoriesJson(counts);
        entry["fraction_sigma"] = categoriesJson(errors.fraction_sigma);
        entry["fraction_sigma_relative"] = categoriesJson(relative);
        entry["hp_np_events"] = errors.hp_np_events;
        entry["non_hp_np_events"] = n * lib.non_hp_np_probability;
        entry["non_hp_np_sigma"] = std::sqrt(n * lib.non_hp_variance);
        entry["non_hp_share"] =
            lib.non_hp_np_probability / (lib.non_hp_np_probability + lib.hp_np_probability);
        entry["legendre_expected"] = lib.legendre;
        entry["legendre_spread"] = lib.legendre_spread;
        entry["legendre_sigma"] = errors.legendre_sigma;
        entry["u_min_expected_count"] =
            errors.hp_np_events *
            *std::min_element(lib.u_probabilities.begin(), lib.u_probabilities.end());
        if (m > 0) {
            const PredictionShift shift = predictionShift(run.libraries.front(), lib,
                                                          expectedErrors(run, 0), run.events);
            entry["shift_from_first_in_sigma"] = {{"fractions", shift.fractions},
                                                  {"legendre", shift.legendre},
                                                  {"u_noncentrality", shift.u_noncentrality}};
            std::vector<double> legendre_shift;
            const ExpectedErrors reference = expectedErrors(run, 0);
            for (std::size_t l = 0; l < lib.legendre.size(); ++l) {
                legendre_shift.push_back((lib.legendre.at(l) - run.libraries.front().legendre.at(l)) /
                                         reference.legendre_sigma.at(l));
            }
            entry["legendre_shift_in_sigma"] = legendre_shift;
            std::array<double, kFirstInteractionCategories> fraction_shift{};
            for (std::size_t c = 0; c < kFirstInteractionCategories; ++c) {
                const double sigma = reference.fraction_sigma.at(c);
                fraction_shift.at(c) =
                    sigma > 0.0 ? (lib.probabilities.at(c) -
                                   run.libraries.front().probabilities.at(c)) /
                                      sigma
                                : 0.0;
            }
            entry["fraction_shift_in_sigma"] = categoriesJson(fraction_shift);
        }
        doc["libraries"].push_back(std::move(entry));
    }
    return doc;
}

} // namespace riptide::sim::checks
