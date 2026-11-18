#pragma once

// Cross-section files of G4NDL, the data library of the Geant4
// high-precision neutron models, read as G4ParticleHPIsoData::Init reads
// them, so that the checks can compare the data Geant4 actually uses with
// the evaluated files it is derived from.
//
// The text (after zlib decompression) is: two integers the model ignores,
// the number of points N, then N pairs (E [eV], sigma [barn]). The model
// stores the points in a G4ParticleHPVector, interpolated lin-lin unless an
// interpolation record is read first, which the cross-section files lack.

#include <istream>
#include <stdexcept>
#include <vector>

namespace riptide::sim::checks {

class HpDataError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Tabulated cross section in riptide::units: En in MeV, sigma in mm^2.
struct HpCrossSection {
    std::vector<double> En;
    std::vector<double> sigma;
};

// Reads one file. An energy may appear twice in a row, marking a
// discontinuity (the two values below and above it). Throws HpDataError if
// the text is truncated or malformed, if N <= 0, if the energies decrease or
// repeat more than twice, or if a cross section is negative.
[[nodiscard]] HpCrossSection parseHpCrossSection(std::istream& in);

// Lin-lin interpolation at En [MeV]; at a discontinuity, the value above
// it. Throws HpDataError outside [En.front(), En.back()]: the checks never
// rely on what the model does beyond its table.
[[nodiscard]] double evaluate(const HpCrossSection& table, double En);

} // namespace riptide::sim::checks
