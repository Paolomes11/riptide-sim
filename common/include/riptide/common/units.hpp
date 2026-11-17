#pragma once

// Unit system of the whole project: MeV, mm, ns, rad.
//
// A quantity is stored as a plain double already expressed in these units; the
// constants below only make conversions explicit at the boundaries (input files,
// printing). Geant4's internal system (CLHEP) uses the same four base units, so
// values cross the sim/analysis boundary unchanged; sim/tests checks this.

namespace riptide::units {

// NOLINTBEGIN(readability-identifier-naming): physics unit symbols keep their usual case
inline constexpr double MeV = 1.0;
inline constexpr double eV = 1.0e-6 * MeV;
inline constexpr double keV = 1.0e-3 * MeV;
inline constexpr double GeV = 1.0e3 * MeV;

inline constexpr double mm = 1.0;
inline constexpr double cm = 10.0 * mm;
inline constexpr double m = 1000.0 * mm;
inline constexpr double mm2 = mm * mm;
inline constexpr double m2 = m * m;
// Evaluated nuclear data give cross sections in barn.
inline constexpr double barn = 1.0e-28 * m2;

inline constexpr double ns = 1.0;

// Written out because M_PI is a POSIX extension, not standard C++17.
inline constexpr double pi = 3.14159265358979323846;
inline constexpr double rad = 1.0;
inline constexpr double deg = pi / 180.0 * rad;
// NOLINTEND(readability-identifier-naming)

} // namespace riptide::units
