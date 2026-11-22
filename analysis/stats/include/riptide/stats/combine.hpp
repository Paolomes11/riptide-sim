#pragma once

// Combination of conditionally independent events sharing one parameter:
// log_post = log_prior + sum_k loglik_k. The prior enters once, whatever the
// number of events; adding it per event would weigh it N times.

#include <vector>

#include "riptide/stats/case_a.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// Throws std::invalid_argument if log_prior does not match the matrix width.
[[nodiscard]] std::vector<double> combineLogLikelihoods(const LogLikelihoodMatrix& loglik,
                                                        const std::vector<double>& log_prior);

} // namespace riptide::stats
