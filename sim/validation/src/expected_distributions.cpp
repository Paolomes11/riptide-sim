#include "riptide/sim/expected_distributions.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "riptide/common/np_kinematics.hpp"

namespace riptide::sim::checks {

namespace {

using schema::FirstInteraction;

// Cross sections read from a file are consistent up to the rounding of their
// six-digit fields; a larger excess of the partial over the total is an error.
constexpr double kPartialSumTolerance = 1e-5;

[[noreturn]] void invalid(const std::string& what) {
    throw std::invalid_argument(what);
}

// Bin edges must be finite, strictly increasing and inside [low, high].
void checkEdges(const std::vector<double>& edges, double low, double high, const char* where) {
    if (edges.size() < 2) {
        invalid(std::string{where} + ": at least two bin edges are needed");
    }
    // Negated range test so that NaN is refused too.
    // NOLINTNEXTLINE(readability-simplify-boolean-expr)
    if (!(edges.front() >= low && edges.back() <= high)) {
        invalid(std::string{where} + ": bin edges outside [" + std::to_string(low) + ", " +
                std::to_string(high) + "]");
    }
    for (std::size_t i = 1; i < edges.size(); ++i) {
        if (!(edges[i] > edges[i - 1])) {
            invalid(std::string{where} + ": bin edges must increase");
        }
    }
}

std::size_t index(FirstInteraction category) {
    return static_cast<std::size_t>(category);
}

} // namespace

double legendreCdf(const std::vector<double>& a, double mu) {
    // NOLINTNEXTLINE(readability-simplify-boolean-expr)
    if (!(mu >= -1.0 && mu <= 1.0)) {
        invalid("legendreCdf: mu = " + std::to_string(mu) + " outside [-1, 1]");
    }
    double cdf = 0.5 * (mu + 1.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        const int degree = static_cast<int>(i) + 1;
        cdf += 0.5 * a[i] *
               (kinematics::legendreP(degree + 1, mu) - kinematics::legendreP(degree - 1, mu));
    }
    return cdf;
}

std::vector<double> protonEnergyBinProbabilities(const std::vector<double>& a,
                                                 const std::vector<double>& edges) {
    checkEdges(edges, 0.0, 1.0, "protonEnergyBinProbabilities");
    // mu = 1 - 2u decreases with u, so the bin [u1, u2] maps to [mu(u2), mu(u1)].
    const auto cosine = [](double u) { return 1.0 - 2.0 * u; };
    std::vector<double> probabilities;
    probabilities.reserve(edges.size() - 1);
    for (std::size_t i = 1; i < edges.size(); ++i) {
        probabilities.push_back(legendreCdf(a, cosine(edges[i - 1])) -
                                legendreCdf(a, cosine(edges[i])));
    }
    return probabilities;
}

std::vector<double> truncatedExponentialBinProbabilities(double sigma_total, double L,
                                                         const std::vector<double>& edges) {
    if (!(sigma_total >= 0.0) || !std::isfinite(sigma_total)) {
        invalid("truncatedExponentialBinProbabilities: Sigma must be finite and >= 0");
    }
    if (!(L > 0.0) || !std::isfinite(L)) {
        invalid("truncatedExponentialBinProbabilities: L must be finite and > 0");
    }
    checkEdges(edges, 0.0, L, "truncatedExponentialBinProbabilities");
    std::vector<double> probabilities;
    probabilities.reserve(edges.size() - 1);
    if (sigma_total <= 0.0) {
        for (std::size_t i = 1; i < edges.size(); ++i) {
            probabilities.push_back((edges[i] - edges[i - 1]) / L);
        }
        return probabilities;
    }
    // e^{-S z1} - e^{-S z2} over 1 - e^{-S L}, written with expm1 so that a
    // thin slab (S L << 1) keeps its digits instead of losing them to 1 - 1.
    const double norm = -std::expm1(-sigma_total * L);
    for (std::size_t i = 1; i < edges.size(); ++i) {
        probabilities.push_back(
            (std::expm1(-sigma_total * edges[i - 1]) - std::expm1(-sigma_total * edges[i])) / norm);
    }
    return probabilities;
}

FirstInteractionProbabilities firstInteraction(const std::vector<MediumComponent>& components,
                                               double L) {
    if (components.empty()) {
        invalid("firstInteraction: the medium has no components");
    }
    if (!(L > 0.0) || !std::isfinite(L)) {
        invalid("firstInteraction: L must be finite and > 0");
    }
    // Macroscopic cross sections [1/mm] of each category.
    std::array<double, kFirstInteractionCategories> partial{};
    double total = 0.0;
    for (const MediumComponent& component : components) {
        const double n = component.number_density;
        const double elastic = component.sigma_elastic;
        const double capture = component.sigma_capture;
        const double sigma = component.sigma_total;
        // Negated test so that NaN is refused too.
        // NOLINTNEXTLINE(readability-simplify-boolean-expr)
        if (!(n >= 0.0 && elastic >= 0.0 && capture >= 0.0 && sigma >= 0.0)) {
            invalid("firstInteraction: densities and cross sections must be >= 0");
        }
        double nonelastic = sigma - elastic - capture;
        if (nonelastic < -kPartialSumTolerance * sigma) {
            invalid("firstInteraction: elastic + capture exceed the total cross section");
        }
        nonelastic = std::max(nonelastic, 0.0);
        const bool hydrogen = component.target == Target::Hydrogen;
        partial.at(index(hydrogen ? FirstInteraction::NpElastic : FirstInteraction::NCElastic)) +=
            n * elastic;
        partial.at(index(hydrogen ? FirstInteraction::Other : FirstInteraction::NCInelastic)) +=
            n * nonelastic;
        partial.at(index(FirstInteraction::Other)) += n * capture;
        total += n * sigma;
    }
    if (!(total > 0.0) || !std::isfinite(total)) {
        invalid("firstInteraction: the total macroscopic cross section must be > 0");
    }
    FirstInteractionProbabilities result;
    result.sigma_total = total;
    // The partial sums are renormalized to their own total, which differs
    // from `total` only by the clamped rounding excess.
    double partial_sum = 0.0;
    for (const double value : partial) {
        partial_sum += value;
    }
    const double interacts = -std::expm1(-total * L);
    for (std::size_t i = 0; i < kFirstInteractionCategories; ++i) {
        result.conditional.at(i) = partial.at(i) / partial_sum;
        result.unconditional.at(i) = interacts * result.conditional.at(i);
    }
    result.unconditional.at(index(FirstInteraction::Unknown)) = std::exp(-total * L);
    return result;
}

} // namespace riptide::sim::checks
