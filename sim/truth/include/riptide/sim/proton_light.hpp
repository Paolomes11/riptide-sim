#pragma once

// Birks light of a proton as a function of its kinetic energy, for the
// visible energy of each proton step in the simulation (D21).
//
// Applying Birks' law step by step, eloss / (1 + kB * eloss / length), makes
// the light of a stopping proton depend on the step function: the quenching
// is concave in dE/dx, so long steps give more light, while short steps with
// fluctuating losses give less (D21). The light of a step is taken instead as
// L(K_pre) - L(K_pre - Edep), with
//     L(K) = integral from 0 to K of dE / (1 + kB * S(E)),
// the integral of the mean stopping power S that drives the transport. The
// light of a proton stopping from Ep is then L(Ep) whatever the steps were,
// which is how kB is defined, while the deposit keeps its fluctuations.
//
// The table is Geant4-free: the simulation fills it with Geant4's stopping
// power at energies of its choice.

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace riptide::sim {

class ProtonLightError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Simpson panels per interval of the stopping-power table. The table of the
// simulation is dense (D21), so the integrand changes by a few percent at
// most across an interval and this many panels leave a relative error far
// below the 1e-6 the tests require.
inline constexpr int kProtonLightPanels = 8;

class ProtonLightTable {
public:
    // K: kinetic energies [MeV], strictly increasing and positive, at least
    // two; S: the stopping power [MeV/mm] at each, finite and positive; kB
    // [mm/MeV] finite and not negative. Between two energies S is log-log;
    // below the first one it is held at its first value. Throws
    // ProtonLightError otherwise.
    ProtonLightTable(std::vector<double> K, std::vector<double> S, double kB);

    // L(K) [MeV] for K in [0, maxEnergy()]; throws ProtonLightError outside.
    [[nodiscard]] double light(double K) const;

    // Visible energy of a step of a proton entering it with K_pre [MeV] and
    // depositing edep [MeV]: L(K_pre) - L(K_pre - edep). A deposit larger
    // than K_pre (the end-of-range overshoot of nuclear stopping, or a local
    // deposit of a hadronic interaction) has its excess quenched as the
    // lowest energies are, so with kB = 0 the result is edep exactly. The
    // result is clamped to [0, edep] against rounding. Throws
    // ProtonLightError if K_pre is outside [0, maxEnergy()] or edep is
    // negative or not finite.
    [[nodiscard]] double visibleEnergy(double K_pre, double edep) const;

    [[nodiscard]] double maxEnergy() const {
        return K_.back();
    }
    [[nodiscard]] double kB() const {
        return kB_;
    }

private:
    [[nodiscard]] double quench(double S) const;
    [[nodiscard]] double intervalLight(std::size_t i, double upper) const;

    std::vector<double> K_;
    std::vector<double> S_;
    double kB_;
    std::vector<double> cumulative_; // L at each tabulated energy
};

} // namespace riptide::sim
