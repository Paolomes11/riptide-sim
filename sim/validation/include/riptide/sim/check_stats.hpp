#pragma once

// Statistics of the physics checks: Pearson chi^2 of a histogram against its
// expected counts, its p-value, the comparison of category fractions with
// their expected probabilities, and the Legendre moments of a sample of
// cosines. Pure functions without Geant4 or ROOT, so they are unit-tested on
// synthetic samples in every preset.

#include <cstddef>
#include <cstdint>
#include <vector>

namespace riptide::sim::checks {

// Smallest expected count accepted in a chi^2 bin or category. Below it the
// Pearson statistic is no longer chi^2-distributed and its p-value would be
// misleading, so the caller must merge bins instead.
inline constexpr double kMinExpectedCount = 5.0;

struct HistogramComparison {
    double chi2 = 0.0;
    int ndf = 0;
    double p_value = 0.0;
    double min_expected = 0.0;
    std::vector<double> pulls; // (o - e) / sqrt(e), one per bin
};

// Pearson chi^2 = sum (o - e)^2 / e of observed against expected counts.
// `constraints` is the number of conditions the expected counts already
// satisfy (1 when they are normalized to the observed total), so
// ndf = bins - constraints. Throws std::invalid_argument on size mismatch,
// negative observed counts, an expected count below kMinExpectedCount, or
// ndf < 1.
[[nodiscard]] HistogramComparison compareHistogram(const std::vector<double>& observed,
                                                   const std::vector<double>& expected,
                                                   int constraints);

// P(X >= chi2) for X ~ chi^2(ndf): the regularized upper incomplete gamma
// function Q(ndf/2, chi2/2). Throws std::invalid_argument for ndf < 1 or a
// negative or non-finite chi2.
[[nodiscard]] double chi2SurvivalProbability(double chi2, int ndf);

struct MergedBins {
    std::vector<double> observed;
    std::vector<double> expected;
    std::vector<double> edges; // observed.size() + 1 bin edges
};

// Merges adjacent bins, left to right, until each merged bin has an expected
// count of at least `min_expected`; a short remainder at the upper end is
// folded into the last complete bin. Deterministic and order-preserving, so
// the same histogram always yields the same chi^2 bins, and the totals are
// conserved exactly. It is applied before compareHistogram, whose p-value is
// meaningless below kMinExpectedCount. Throws std::invalid_argument on size
// mismatch (edges must have one entry more than the bins), no bins, a
// non-positive min_expected, or a total expected count below min_expected.
[[nodiscard]] MergedBins mergeSparseBins(const std::vector<double>& observed,
                                         const std::vector<double>& expected,
                                         const std::vector<double>& edges, double min_expected);

struct FractionComparison {
    std::uint64_t total = 0;
    double chi2 = 0.0;
    int ndf = 0;
    double p_value = 0.0;
    std::vector<double> fractions; // k_i / n
    std::vector<double> sigmas;    // sqrt(p_i (1 - p_i) / n), from the expected p_i
    std::vector<double> pulls;     // (k_i / n - p_i) / sigma_i
};

// Counts k_i of mutually exclusive categories against their expected
// probabilities p_i (summing to 1). The chi^2 is the multinomial Pearson
// statistic with ndf = categories - 1; the per-category pulls are binomial
// and are not independent, so their squares do not add up to chi2. The
// sigma uses the expected p_i, which keeps it meaningful when k_i = 0.
// Throws std::invalid_argument on size mismatch, fewer than two categories,
// probabilities outside (0, 1) or not summing to 1, or an expected count
// n p_i below kMinExpectedCount.
[[nodiscard]] FractionComparison compareFractions(const std::vector<std::uint64_t>& counts,
                                                  const std::vector<double>& probabilities);

struct LegendreMoments {
    std::size_t samples = 0;
    std::vector<double> values; // values[l-1] = <P_l(mu)>, l = 1 ... max_degree
    std::vector<double> errors; // standard error of each mean
};

// Sample estimate of the Legendre coefficients of an angular distribution.
// For f(mu) = sum_l (2l+1)/2 a_l P_l(mu), orthogonality gives
// <P_l(mu)> = a_l, so the sample mean of P_l estimates a_l with standard
// error s_l / sqrt(N), s_l^2 being the unbiased sample variance of P_l.
// Throws std::invalid_argument for max_degree < 1, fewer than two samples,
// or a cosine outside [-1, 1] (a sign of a bug upstream, not to be clamped).
[[nodiscard]] LegendreMoments legendreMoments(const std::vector<double>& mu, int max_degree);

} // namespace riptide::sim::checks
