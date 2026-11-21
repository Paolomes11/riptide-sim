#include "riptide/sim/pstar.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "riptide/common/units.hpp"

namespace riptide::sim::checks {

namespace {

constexpr std::size_t kColumns = 7;

// The table gives four significant digits, so each of the three stopping
// powers carries a relative rounding error up to 5e-4 of its own value: the
// sum check allows three such errors on the total.
constexpr double kSumTolerance = 1.5e-3;

// Mass units of the table, expressed in riptide::units for a density in
// g/cm3: (MeV cm2/g) * (g/cm3) = MeV/cm, and (g/cm2) / (g/cm3) = cm.
constexpr double kStoppingUnit = units::MeV / units::cm;
constexpr double kRangeUnit = units::cm;

// Parses a line as exactly seven numbers; false for anything else.
bool parseRow(const std::string& line, std::array<double, kColumns>& row) {
    std::istringstream fields(line);
    for (double& value : row) {
        if (!(fields >> value)) {
            return false;
        }
    }
    std::string extra;
    return !(fields >> extra);
}

bool isBlank(const std::string& line) {
    return line.find_first_not_of(" \t\r") == std::string::npos;
}

// Index lo of the interval [Ep[lo], Ep[lo+1]] holding Ep, for Ep inside
// the table and strictly below its last energy.
std::size_t intervalBelow(const std::vector<double>& energies, double Ep) {
    auto upper = std::upper_bound(energies.begin(), energies.end(), Ep);
    return static_cast<std::size_t>(std::distance(energies.begin(), upper)) - 1;
}

// Power law through (x0, y0) and (x1, y1), all positive, evaluated at x.
double logLog(double x0, double y0, double x1, double y1, double x) {
    const double slope = std::log(y1 / y0) / std::log(x1 / x0);
    return y0 * std::exp(slope * std::log(x / x0));
}

double interpolate(const std::vector<double>& energies, const std::vector<double>& values,
                   double Ep, const char* what) {
    if (energies.empty() || !(Ep >= energies.front()) || !(Ep <= energies.back())) {
        throw PstarError(std::string("PSTAR ") + what + ": energy " + std::to_string(Ep) +
                         " MeV outside the table");
    }
    if (!(Ep < energies.back())) {
        return values.back();
    }
    const std::size_t lo = intervalBelow(energies, Ep);
    return logLog(energies[lo], values[lo], energies[lo + 1], values[lo + 1], Ep);
}

} // namespace

PstarTable parsePstarTable(std::istream& in, double density_g_per_cm3) {
    if (!(density_g_per_cm3 > 0.0) || !std::isfinite(density_g_per_cm3)) {
        throw PstarError("PSTAR: density must be finite and positive");
    }
    PstarTable table;
    bool seen_title = false;
    bool in_data = false;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        std::array<double, kColumns> row{};
        if (!parseRow(line, row)) {
            if (in_data && !isBlank(line)) {
                throw PstarError("PSTAR: malformed data row at line " +
                                 std::to_string(line_number));
            }
            seen_title = seen_title || line.find("PSTAR") != std::string::npos;
            continue;
        }
        if (!seen_title) {
            throw PstarError("PSTAR: data before a header naming PSTAR");
        }
        in_data = true;
        const bool valid = std::all_of(row.begin(), row.end(),
                                       [](double v) { return std::isfinite(v) && v > 0.0; });
        if (!valid) {
            throw PstarError("PSTAR: values must be finite and positive at line " +
                             std::to_string(line_number));
        }
        const auto [energy, s_el, s_nuc, s_tot, csda, projected, detour] = row;
        static_cast<void>(detour);
        if (!table.Ep.empty() && !(energy * units::MeV > table.Ep.back())) {
            throw PstarError("PSTAR: energies not increasing at line " +
                             std::to_string(line_number));
        }
        if (std::abs(s_tot - (s_el + s_nuc)) > kSumTolerance * s_tot) {
            throw PstarError("PSTAR: total stopping power is not the sum of its parts at line " +
                             std::to_string(line_number));
        }
        const double stopping = density_g_per_cm3 * kStoppingUnit;
        const double range = kRangeUnit / density_g_per_cm3;
        table.Ep.push_back(energy * units::MeV);
        table.S_el.push_back(s_el * stopping);
        table.S_nuc.push_back(s_nuc * stopping);
        table.S_tot.push_back(s_tot * stopping);
        table.csda_range.push_back(csda * range);
        table.projected_range.push_back(projected * range);
    }
    if (table.Ep.empty()) {
        throw PstarError("PSTAR: no data rows");
    }
    return table;
}

double stoppingPower(const PstarTable& table, double Ep) {
    return interpolate(table.Ep, table.S_tot, Ep, "stopping power");
}

double csdaRange(const PstarTable& table, double Ep) {
    return interpolate(table.Ep, table.csda_range, Ep, "CSDA range");
}

namespace {

void checkLightArguments(const PstarTable& table, double kB, int panels) {
    if (table.Ep.empty()) {
        throw PstarError("PSTAR light: empty table");
    }
    if (!(kB >= 0.0) || !std::isfinite(kB)) {
        throw PstarError("PSTAR light: kB must be finite and not negative");
    }
    if (panels < 1) {
        throw PstarError("PSTAR light: at least one Simpson panel per interval");
    }
}

void checkLightEnergy(const PstarTable& table, double Ep) {
    if (!(Ep >= 0.0) || !(Ep <= table.Ep.back())) {
        throw PstarError("PSTAR light: energy " + std::to_string(Ep) +
                         " MeV negative or above the table");
    }
}

double quench(double kB, double S) {
    return 1.0 / (1.0 + kB * S);
}

// Light below the first tabulated energy, where the stopping power is
// unknown: S is held at its first value, an error bounded by that energy.
double lightBelowTable(const PstarTable& table, double Ep, double kB) {
    return Ep * quench(kB, table.S_tot.front());
}

// Simpson's rule from Ep[i] to `upper` <= Ep[i+1], with S log-log between
// the ends of the interval, so the integrand is smooth inside every panel.
double intervalLight(const PstarTable& table, std::size_t i, double upper, double kB, int panels) {
    const double e0 = table.Ep[i];
    const double e1 = table.Ep[i + 1];
    const double s0 = table.S_tot[i];
    const double s1 = table.S_tot[i + 1];
    const auto integrand = [&](double E) { return quench(kB, logLog(e0, s0, e1, s1, E)); };
    const double width = (upper - e0) / panels;
    double sum = 0.0;
    for (int k = 0; k < panels; ++k) {
        const double a = e0 + k * width;
        const double b = (k + 1 == panels) ? upper : a + width;
        sum += (b - a) * (integrand(a) + 4.0 * integrand(0.5 * (a + b)) + integrand(b));
    }
    return sum / 6.0;
}

} // namespace

double visibleEnergy(const PstarTable& table, double Ep, double kB, int panels) {
    checkLightArguments(table, kB, panels);
    checkLightEnergy(table, Ep);
    const double first = table.Ep.front();
    if (!(Ep > first)) {
        return lightBelowTable(table, Ep, kB);
    }
    double light = lightBelowTable(table, first, kB);
    for (std::size_t i = 0; i + 1 < table.Ep.size() && table.Ep[i] < Ep; ++i) {
        light += intervalLight(table, i, std::min(table.Ep[i + 1], Ep), kB, panels);
    }
    return light;
}

BirksLightTable::BirksLightTable(PstarTable table, double kB, int panels)
    : table_(std::move(table)), kB_(kB), panels_(panels) {
    checkLightArguments(table_, kB_, panels_);
    cumulative_.reserve(table_.Ep.size());
    cumulative_.push_back(lightBelowTable(table_, table_.Ep.front(), kB_));
    for (std::size_t i = 0; i + 1 < table_.Ep.size(); ++i) {
        cumulative_.push_back(cumulative_.back() +
                              intervalLight(table_, i, table_.Ep[i + 1], kB_, panels_));
    }
}

double BirksLightTable::light(double Ep) const {
    checkLightEnergy(table_, Ep);
    if (!(Ep > table_.Ep.front())) {
        return lightBelowTable(table_, Ep, kB_);
    }
    // Ep in (Ep[lo], Ep[lo+1]]: the full intervals below lo are in the sum.
    const auto upper = std::lower_bound(table_.Ep.begin(), table_.Ep.end(), Ep);
    const auto lo = static_cast<std::size_t>(std::distance(table_.Ep.begin(), upper)) - 1;
    return cumulative_[lo] + intervalLight(table_, lo, Ep, kB_, panels_);
}

} // namespace riptide::sim::checks
