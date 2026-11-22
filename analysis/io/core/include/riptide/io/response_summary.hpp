#pragma once

// Summaries of the detector response applied to a sample, written to
// analysis.json: they show that the HP conversion and the smearing did what
// they claim on the events actually used, not only on the unit-test fixtures.
//
// Smearing residuals are judged against the known widths:
//  - Ep and theta_p residuals are N(0, sigma): the mean has the pull
//    mean / (sigma / sqrt(N)), and the sample variance s^2 (N - 1 in the
//    denominator) the pull (s^2 / sigma^2 - 1) / sqrt(2 / (N - 1)), since
//    (N - 1) s^2 / sigma^2 is chi^2 with N - 1 degrees of freedom;
//  - the angle delta between the true and the smeared track is the norm of a
//    2D Gaussian offset, so delta^2 / (2 sigma^2) is Exp(1), with mean and
//    standard deviation 1: its mean has the pull (mean - 1) sqrt(N).
// With a zero width the residuals must vanish; only their largest absolute
// value is reported and the pulls are null.

#include <nlohmann/json.hpp>

#include "riptide/common/event_sample.hpp"
#include "riptide/io/response.hpp"

namespace riptide::io {

// Relative offset (truth_Ep - truth_Ep_phys) / truth_Ep_phys: min, max, mean;
// and the largest relative difference between truth_Ep_phys and the
// bare-proton energy at the true lab angle. Needs truth_Ep_phys filled.
// Events with truth_Ep_phys = 0 (a recoil at exactly 90 degrees) have no
// relative offset and are counted apart. Throws IoError on an empty sample or
// missing columns.
[[nodiscard]] nlohmann::json hpConversionSummary(const EventSample& sample);

// Residual summaries of reco_Ep, reco_theta_p and reco_p_dir against the
// truth they were drawn from. Needs at least two events and every reco
// column filled. Throws IoError otherwise.
[[nodiscard]] nlohmann::json resolutionSummary(const EventSample& sample,
                                               const Resolution& resolution);

} // namespace riptide::io
