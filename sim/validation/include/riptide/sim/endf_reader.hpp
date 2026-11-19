#pragma once

// Reader for the parts of an ENDF-6 neutron file (ENDF-102 manual) that the
// physics checks need: cross sections (MF3) and Legendre angular
// distributions (MF4, LTT=1). Values are converted to the units of
// riptide::units at reading time (MeV, mm^2), so no caller ever sees eV or barn.

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <istream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "riptide/sim/endf_format.hpp"

namespace riptide::sim::endf {

// The records of one (MF, MT) section, without its SEND record.
struct EndfSection {
    std::size_t first_line = 0; // 1-based line number in the source, for messages
    std::vector<EndfLine> lines;
};

// One material of an ENDF tape. Only single-material files are accepted, as
// distributed by the evaluation projects for each nuclide.
struct EndfMaterial {
    std::string source; // file name or caller label, for messages
    int mat = 0;
    double za = 0.0;  // 1000 Z + A, from the MF1/MT451 HEAD record
    double awr = 0.0; // mass in units of the neutron mass
    std::map<std::pair<int, int>, EndfSection> sections; // key (MF, MT)
};

// Groups the lines of a tape by (MF, MT). The tape identification line and
// the SEND/FEND/MEND/TEND end records are dropped. Throws EndfError (with
// source and line number) on malformed lines, a second material, a section
// split in two, or a missing MF1/MT451.
[[nodiscard]] EndfMaterial readEndf(std::istream& input, const std::string& source);
[[nodiscard]] EndfMaterial readEndfFile(const std::filesystem::path& path);

// Cross section of reaction `mt` from MF3, x in MeV and y in mm^2.
[[nodiscard]] Tab1 crossSection(const EndfMaterial& material, int mt);

// Reference frame of an angular distribution (LCT of MF4).
enum class Frame : std::uint8_t {
    Lab = 1,
    CenterOfMass = 2,
};

// Legendre coefficients a_l(E), l >= 1, of the normalized angular
// distribution f(mu, E) = sum_l (2l+1)/2 a_l(E) P_l(mu), with a_0 = 1 implicit.
// The incident-energy grid follows the TAB2 record of the file.
struct LegendreTable {
    Frame frame = Frame::CenterOfMass;
    std::vector<int> boundaries;
    std::vector<InterpolationLaw> laws;
    std::vector<double> energies;                  // MeV
    std::vector<std::vector<double>> coefficients; // coefficients[i][l-1] = a_l(E_i)
};

// MF4 section of reaction `mt` (default: elastic scattering). Only LTT=1
// (Legendre coefficients over the whole energy range) with LI=0 is accepted;
// the other representations throw EndfError.
[[nodiscard]] LegendreTable legendreAngular(const EndfMaterial& material, int mt = 2);

// a_1 ... a_L at incident energy En [MeV], interpolated lin-lin between the
// two neighbouring energies. When the number of coefficients differs between
// them, the missing ones are zero, so L is the larger of the two. Throws
// EndfError outside the tabulated energy range.
[[nodiscard]] std::vector<double> coefficients(const LegendreTable& table, double En);

} // namespace riptide::sim::endf
