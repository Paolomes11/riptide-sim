#pragma once

// NIST PSTAR proton stopping-power and range table, the reference for the
// range and the Birks light of the recoil protons in the checks.
//
// The text table has a header (title, material, column names and units)
// followed by rows of seven numbers: E [MeV], electronic, nuclear and total
// stopping power [MeV cm2/g], CSDA and projected range [g/cm2], detour
// factor. The reader converts mass units to linear ones with the density
// the caller gives, so the table is in riptide::units like the simulation.

#include <istream>
#include <stdexcept>
#include <vector>

namespace riptide::sim::checks {

class PstarError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Density of the PSTAR material 216 "PLASTIC SCINTILLATOR
// (vinyltoluene-based)", equal to that of G4_PLASTIC_SC_VINYLTOLUENE in
// G4NistMaterialBuilder, in g/cm3: grams are not a riptide unit.
inline constexpr double kVinyltolueneDensityGramPerCm3 = 1.032;

// One entry per tabulated energy: Ep in MeV, stopping powers in MeV/mm,
// ranges in mm. The detour factor is not kept: it is projected / CSDA.
struct PstarTable {
    std::vector<double> Ep;
    std::vector<double> S_el;
    std::vector<double> S_nuc;
    std::vector<double> S_tot;
    std::vector<double> csda_range;
    std::vector<double> projected_range;
};

// Reads the text table. Lines before the first row of seven numbers are the
// header, which must mention "PSTAR"; after it only data rows or blank lines.
// Throws PstarError on a malformed row, an empty table, energies that do not
// increase, values that are not finite and positive, a non-positive density,
// or a total stopping power that is not the sum of the two parts within the
// rounding of four significant digits.
[[nodiscard]] PstarTable parsePstarTable(std::istream& in, double density_g_per_cm3);

// Total stopping power [MeV/mm] at Ep [MeV], log-log interpolation. Throws
// PstarError outside [Ep.front(), Ep.back()].
[[nodiscard]] double stoppingPower(const PstarTable& table, double Ep);

// CSDA range [mm] at Ep [MeV], log-log interpolation. Throws PstarError
// outside [Ep.front(), Ep.back()].
[[nodiscard]] double csdaRange(const PstarTable& table, double Ep);

// Number of Simpson panels per table interval in visibleEnergy: on the NIST
// table doubling it changes L by less than 1e-9 in relative terms (checked
// by the reference test), far below any statistical error of the checks.
inline constexpr int kBirksPanels = 16;

// Light expected from Birks' law for a proton stopping from Ep [MeV]:
// L = integral from 0 to Ep of dE / (1 + kB * S_tot(E)), kB in mm/MeV, in
// MeV of visible energy. Below the first tabulated energy S is taken
// constant at its first value; the error this makes is below that energy
// (1 keV in the NIST table). Throws PstarError if Ep is above the table, if
// Ep or kB is negative or not finite, or if panels < 1.
[[nodiscard]] double visibleEnergy(const PstarTable& table, double Ep, double kB,
                                   int panels = kBirksPanels);

} // namespace riptide::sim::checks
