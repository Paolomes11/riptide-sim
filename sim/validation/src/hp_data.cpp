#include "riptide/sim/hp_data.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <string>

#include "riptide/common/units.hpp"

namespace riptide::sim::checks {

HpCrossSection parseHpCrossSection(std::istream& in) {
    int ignored_first = 0;
    int ignored_second = 0;
    int n_points = 0;
    if (!(in >> ignored_first >> ignored_second >> n_points)) {
        throw HpDataError("HP cross section: expected two integers and the number of points");
    }
    if (n_points <= 0) {
        throw HpDataError("HP cross section: number of points must be positive, got " +
                          std::to_string(n_points));
    }
    HpCrossSection table;
    table.En.reserve(static_cast<std::size_t>(n_points));
    table.sigma.reserve(static_cast<std::size_t>(n_points));
    for (int i = 0; i < n_points; ++i) {
        double energy_eV = 0.0;
        double sigma_barn = 0.0;
        if (!(in >> energy_eV >> sigma_barn)) {
            throw HpDataError("HP cross section: truncated at point " + std::to_string(i) + " of " +
                              std::to_string(n_points));
        }
        if (!std::isfinite(energy_eV) || !std::isfinite(sigma_barn) || sigma_barn < 0.0) {
            throw HpDataError("HP cross section: invalid value at point " + std::to_string(i));
        }
        const double En = energy_eV * units::eV;
        // A repeated energy marks a discontinuity, as in ENDF tables; a
        // third equal energy or a decrease has no meaning.
        const std::size_t size = table.En.size();
        const bool decreasing = size >= 1 && En < table.En[size - 1];
        const bool third_equal = size >= 2 && !(En > table.En[size - 2]);
        if (decreasing || third_equal) {
            throw HpDataError("HP cross section: energies not increasing at point " +
                              std::to_string(i));
        }
        table.En.push_back(En);
        table.sigma.push_back(sigma_barn * units::barn);
    }
    return table;
}

double evaluate(const HpCrossSection& table, double En) {
    if (table.En.empty() || En < table.En.front() || En > table.En.back()) {
        throw HpDataError("HP cross section: energy " + std::to_string(En) +
                          " MeV outside the table");
    }
    // First point strictly above En, so at a discontinuity the value above
    // it is taken, as endf::evaluate does; at the last point its value.
    auto upper = std::upper_bound(table.En.begin(), table.En.end(), En);
    if (upper == table.En.end()) {
        return table.sigma.back();
    }
    const auto hi = static_cast<std::size_t>(std::distance(table.En.begin(), upper));
    const std::size_t lo = hi - 1;
    const double t = (En - table.En[lo]) / (table.En[hi] - table.En[lo]);
    return table.sigma[lo] + t * (table.sigma[hi] - table.sigma[lo]);
}

} // namespace riptide::sim::checks
