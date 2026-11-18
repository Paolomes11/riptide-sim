#include "riptide/sim/check_stats.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

#include "riptide/common/np_kinematics.hpp"

namespace riptide::sim::checks {

namespace {

// Relative accuracy of the incomplete gamma evaluation; a few units of the
// double epsilon, far below any statistical precision a check can reach.
constexpr double kGammaEpsilon = 4.0 * std::numeric_limits<double>::epsilon();
constexpr int kGammaMaxIterations = 10000;
// Guards the modified Lentz algorithm against a zero denominator.
constexpr double kLentzTiny = std::numeric_limits<double>::min() / kGammaEpsilon;
// Probabilities computed from cross sections sum to 1 up to rounding.
constexpr double kProbabilitySumTolerance = 1e-9;

[[noreturn]] void invalid(const std::string& what) {
    throw std::invalid_argument(what);
}

// P(a, x) by its power series, used for x < a + 1 where it converges fast:
//   P(a, x) = e^{-x} x^a / Gamma(a + 1) * sum_n x^n / ((a + 1) ... (a + n)).
double lowerGammaSeries(double a, double x) {
    double term = 1.0 / a;
    double sum = term;
    for (int n = 1; n <= kGammaMaxIterations; ++n) {
        term *= x / (a + n);
        sum += term;
        if (std::abs(term) < std::abs(sum) * kGammaEpsilon) {
            return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
        }
    }
    throw std::runtime_error("incomplete gamma series did not converge");
}

// Q(a, x) by its continued fraction (modified Lentz), used for x >= a + 1:
//   Q(a, x) = e^{-x} x^a / Gamma(a) * 1/(x+1-a- 1(1-a)/(x+3-a- 2(2-a)/(x+5-a- ...))).
double upperGammaFraction(double a, double x) {
    double b = x + 1.0 - a;
    double c = 1.0 / kLentzTiny;
    double d = 1.0 / b;
    double h = d;
    for (int n = 1; n <= kGammaMaxIterations; ++n) {
        const double an = -n * (n - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < kLentzTiny) {
            d = kLentzTiny;
        }
        c = b + an / c;
        if (std::abs(c) < kLentzTiny) {
            c = kLentzTiny;
        }
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::abs(delta - 1.0) < kGammaEpsilon) {
            return h * std::exp(-x + a * std::log(x) - std::lgamma(a));
        }
    }
    throw std::runtime_error("incomplete gamma continued fraction did not converge");
}

} // namespace

HistogramComparison compareHistogram(const std::vector<double>& observed,
                                     const std::vector<double>& expected, int constraints) {
    if (observed.size() != expected.size()) {
        invalid("compareHistogram: observed and expected differ in size");
    }
    if (constraints < 0) {
        invalid("compareHistogram: negative number of constraints");
    }
    const int ndf = static_cast<int>(observed.size()) - constraints;
    if (ndf < 1) {
        invalid("compareHistogram: ndf = " + std::to_string(ndf) + " < 1");
    }
    HistogramComparison result;
    result.ndf = ndf;
    result.min_expected = std::numeric_limits<double>::infinity();
    result.pulls.reserve(observed.size());
    for (std::size_t i = 0; i < observed.size(); ++i) {
        // Negated tests so that NaN is refused too.
        if (!(observed[i] >= 0.0)) {
            invalid("compareHistogram: bin " + std::to_string(i) + " has a negative count");
        }
        if (!(expected[i] >= kMinExpectedCount)) {
            invalid("compareHistogram: bin " + std::to_string(i) + " expects " +
                    std::to_string(expected[i]) + " < " + std::to_string(kMinExpectedCount) +
                    " entries; merge bins");
        }
        const double pull = (observed[i] - expected[i]) / std::sqrt(expected[i]);
        result.pulls.push_back(pull);
        result.chi2 += pull * pull;
        result.min_expected = std::min(result.min_expected, expected[i]);
    }
    result.p_value = chi2SurvivalProbability(result.chi2, ndf);
    return result;
}

double chi2SurvivalProbability(double chi2, int ndf) {
    if (ndf < 1) {
        invalid("chi2SurvivalProbability: ndf = " + std::to_string(ndf) + " < 1");
    }
    if (!(chi2 >= 0.0) || !std::isfinite(chi2)) {
        invalid("chi2SurvivalProbability: chi2 must be finite and >= 0");
    }
    const double a = 0.5 * ndf;
    const double x = 0.5 * chi2;
    if (x <= 0.0) {
        return 1.0;
    }
    // Each expansion is used where it converges quickly; below a + 1 the
    // series gives P, which is not close to 1 there, so 1 - P keeps its
    // relative accuracy.
    if (x < a + 1.0) {
        return 1.0 - lowerGammaSeries(a, x);
    }
    return upperGammaFraction(a, x);
}

MergedBins mergeSparseBins(const std::vector<double>& observed, const std::vector<double>& expected,
                           const std::vector<double>& edges, double min_expected) {
    if (observed.size() != expected.size() || edges.size() != observed.size() + 1) {
        invalid("mergeSparseBins: observed, expected and edges sizes do not match");
    }
    if (observed.empty()) {
        invalid("mergeSparseBins: no bins");
    }
    if (!(min_expected > 0.0)) {
        invalid("mergeSparseBins: min_expected must be positive");
    }
    MergedBins merged;
    merged.edges.push_back(edges.front());
    double group_observed = 0.0;
    double group_expected = 0.0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        group_observed += observed[i];
        group_expected += expected[i];
        if (group_expected >= min_expected) {
            merged.observed.push_back(group_observed);
            merged.expected.push_back(group_expected);
            merged.edges.push_back(edges[i + 1]);
            group_observed = 0.0;
            group_expected = 0.0;
        }
    }
    if (merged.observed.empty()) {
        invalid("mergeSparseBins: total expected count below min_expected");
    }
    // Leftover upper bins that never reached the threshold join the last group.
    if (group_observed > 0.0 || group_expected > 0.0) {
        merged.observed.back() += group_observed;
        merged.expected.back() += group_expected;
    }
    merged.edges.back() = edges.back();
    return merged;
}

FractionComparison compareFractions(const std::vector<std::uint64_t>& counts,
                                    const std::vector<double>& probabilities) {
    if (counts.size() != probabilities.size()) {
        invalid("compareFractions: counts and probabilities differ in size");
    }
    if (counts.size() < 2) {
        invalid("compareFractions: at least two categories are needed");
    }
    const double sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    if (!(std::abs(sum - 1.0) <= kProbabilitySumTolerance)) {
        invalid("compareFractions: probabilities sum to " + std::to_string(sum) + ", not 1");
    }
    FractionComparison result;
    result.total = std::accumulate(counts.begin(), counts.end(), std::uint64_t{0});
    result.ndf = static_cast<int>(counts.size()) - 1;
    const auto n = static_cast<double>(result.total);
    for (std::size_t i = 0; i < counts.size(); ++i) {
        const double p = probabilities[i];
        // Negated range test so that NaN is refused too.
        // NOLINTNEXTLINE(readability-simplify-boolean-expr)
        if (!(p > 0.0 && p < 1.0)) {
            invalid("compareFractions: probability of category " + std::to_string(i) +
                    " is outside (0, 1)");
        }
        const double expected = n * p;
        if (!(expected >= kMinExpectedCount)) {
            invalid("compareFractions: category " + std::to_string(i) + " expects " +
                    std::to_string(expected) + " < " + std::to_string(kMinExpectedCount) +
                    " entries");
        }
        const auto k = static_cast<double>(counts[i]);
        const double fraction = k / n;
        const double sigma = std::sqrt(p * (1.0 - p) / n);
        result.fractions.push_back(fraction);
        result.sigmas.push_back(sigma);
        result.pulls.push_back((fraction - p) / sigma);
        result.chi2 += (k - expected) * (k - expected) / expected;
    }
    result.p_value = chi2SurvivalProbability(result.chi2, result.ndf);
    return result;
}

LegendreMoments legendreMoments(const std::vector<double>& mu, int max_degree) {
    if (max_degree < 1) {
        invalid("legendreMoments: max_degree must be >= 1");
    }
    if (mu.size() < 2) {
        invalid("legendreMoments: at least two samples are needed");
    }
    for (const double cosine : mu) {
        // NOLINTNEXTLINE(readability-simplify-boolean-expr)
        if (!(cosine >= -1.0 && cosine <= 1.0)) {
            invalid("legendreMoments: cosine " + std::to_string(cosine) + " outside [-1, 1]");
        }
    }
    const auto n = static_cast<double>(mu.size());
    LegendreMoments result;
    result.samples = mu.size();
    for (int degree = 1; degree <= max_degree; ++degree) {
        // Two passes: the mean first, then the squared deviations, which
        // avoids the cancellation of <P^2> - <P>^2 when a_l^2 ~ <P_l^2>.
        double sum = 0.0;
        for (const double cosine : mu) {
            sum += kinematics::legendreP(degree, cosine);
        }
        const double mean = sum / n;
        double squares = 0.0;
        for (const double cosine : mu) {
            const double deviation = kinematics::legendreP(degree, cosine) - mean;
            squares += deviation * deviation;
        }
        result.values.push_back(mean);
        result.errors.push_back(std::sqrt(squares / (n - 1.0) / n));
    }
    return result;
}

} // namespace riptide::sim::checks
