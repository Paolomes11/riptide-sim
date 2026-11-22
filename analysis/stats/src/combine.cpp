#include "riptide/stats/combine.hpp"

#include <stdexcept>

namespace riptide::stats {

std::vector<double> combineLogLikelihoods(const LogLikelihoodMatrix& loglik,
                                          const std::vector<double>& log_prior) {
    if (log_prior.size() != loglik.points) {
        throw std::invalid_argument("the prior does not match the log-likelihood grid");
    }
    std::vector<double> result = log_prior;
    for (std::size_t k = 0; k < loglik.events; ++k) {
        const double* row = loglik.row(k);
        for (std::size_t i = 0; i < loglik.points; ++i) {
            result[i] += row[i];
        }
    }
    return result;
}

} // namespace riptide::stats
