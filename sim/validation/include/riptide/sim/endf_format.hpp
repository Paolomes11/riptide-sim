#pragma once

// Low-level pieces of the ENDF-6 format (ENDF-102 manual, section 0) used to
// read evaluated neutron data: the 80-column record, its six 11-character
// fields, the Fortran number notation, and the one-dimensional interpolation
// laws of TAB1 records. Geant4-free, so it is unit-tested in every preset.

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace riptide::sim::endf {

class EndfError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

inline constexpr std::size_t kFieldWidth = 11;
inline constexpr std::size_t kFieldsPerLine = 6;
// Columns 1-66 hold the data fields, 67-70 MAT, 71-72 MF, 73-75 MT; columns
// 76-80 (line sequence number) are optional and ignored.
inline constexpr std::size_t kDataWidth = kFieldWidth * kFieldsPerLine;
inline constexpr std::size_t kMinLineWidth = 75;

// Parses one data field. Accepts the Fortran forms found in ENDF files, where
// the exponent letter is usually omitted ("9.991673-1", "-1.2+3"), as well as
// "1.0E+00" and plain integers. A blank field is 0, as in Fortran list input.
// Anything else throws EndfError: a silently misread number would corrupt a
// cross section without any visible sign.
[[nodiscard]] double parseEndfFloat(std::string_view field);

// Parses an integer field (right-justified, blank = 0); throws on non-integers.
[[nodiscard]] int parseEndfInt(std::string_view field);

// One record of an ENDF file.
struct EndfLine {
    int mat = 0;
    int mf = 0;
    int mt = 0;
    std::array<std::string, kFieldsPerLine> fields;

    [[nodiscard]] double real(std::size_t index) const {
        return parseEndfFloat(fields.at(index));
    }
    [[nodiscard]] int integer(std::size_t index) const {
        return parseEndfInt(fields.at(index));
    }
};

// Splits one line into fields and MAT/MF/MT; throws if it is shorter than 75
// columns or MAT/MF/MT are not integers.
[[nodiscard]] EndfLine parseEndfLine(std::string_view line);

// Interpolation laws of Table 16 of ENDF-102. Law 6 (charged-particle
// Coulomb penetrability) is not needed for neutron data and is refused.
enum class InterpolationLaw : std::uint8_t {
    Histogram = 1, // y constant, equal to the value at the lower end
    LinLin = 2,    // y linear in x
    LinLog = 3,    // y linear in ln(x)
    LogLin = 4,    // ln(y) linear in x
    LogLog = 5,    // ln(y) linear in ln(x)
};

// Converts the INT code of a file; throws for codes outside 1-5.
[[nodiscard]] InterpolationLaw interpolationLaw(int code);

// Interpolates between (x1, y1) and (x2, y2), x1 < x2, x1 <= x <= x2.
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters): the two points are the physics signature
[[nodiscard]] double interpolate(InterpolationLaw law, double x1, double y1, double x2, double y2,
                                 double x);

// Tabulated function of a TAB1 record. Range m (0-based) covers the points
// with 1-based index up to boundaries[m] (the NBT values of the file) and
// uses laws[m] between consecutive points.
struct Tab1 {
    std::vector<int> boundaries;
    std::vector<InterpolationLaw> laws;
    std::vector<double> x;
    std::vector<double> y;
};

// Checks the invariants of a Tab1 (sizes, NBT strictly increasing and ending at
// the last point, x non-decreasing, at most two equal x in a row); throws
// EndfError with `where` in the message otherwise.
void validate(const Tab1& table, std::string_view where);

// Value of the tabulated function at x. At a discontinuity (x repeated, two
// y values) the value above it is returned, since an interval of zero width
// is skipped. Throws EndfError outside [x.front(), x.back()]: extrapolating
// evaluated data is never meaningful for these checks.
[[nodiscard]] double evaluate(const Tab1& table, double x);

} // namespace riptide::sim::endf
