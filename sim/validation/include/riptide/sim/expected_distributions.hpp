#pragma once

// Expected distributions of the physics checks, computed from evaluated data
// (Legendre coefficients, cross sections) and from the medium. The caller
// passes the numbers in; nothing here reads files, Geant4 or ROOT.

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "riptide/common/event_schema.hpp"

namespace riptide::sim::checks {

// Cumulative distribution of the normalized angular distribution
// f(mu) = 1/2 + sum_l (2l+1)/2 a_l P_l(mu), with a[l-1] = a_l:
//   F(mu) = (mu + 1)/2 + sum_l a_l/2 (P_{l+1}(mu) - P_{l-1}(mu)),
// from the integral of P_l over [-1, mu], (P_{l+1} - P_{l-1})/(2l+1).
// Throws std::invalid_argument for mu outside [-1, 1].
[[nodiscard]] double legendreCdf(const std::vector<double>& a, double mu);

// Probability of each bin of u = Ep / Ep_max for elastic n-p scattering whose
// neutron CM angular distribution has coefficients `a`. The map
// u = (1 - mu) / 2 is exact for elastic scattering (cosCmFromProtonEnergy),
// so a bin [u1, u2] is the cosine interval [1 - 2 u2, 1 - 2 u1]. The edges must
// increase within [0, 1]; the probabilities sum to 1 when they span it.
// Throws std::invalid_argument otherwise.
[[nodiscard]] std::vector<double> protonEnergyBinProbabilities(const std::vector<double>& a,
                                                               const std::vector<double>& edges);

// Probability of each depth bin of the first interaction point, given that an
// interaction happens within the thickness L of a slab crossed at normal
// incidence:
//   p(z) = Sigma e^{-Sigma z} / (1 - e^{-Sigma L}),  0 <= z <= L,
// with Sigma the total macroscopic cross section [1/mm]. Sigma = 0 gives the
// uniform limit. The edges [mm] must increase within [0, L]. Throws
// std::invalid_argument otherwise.
[[nodiscard]] std::vector<double>
truncatedExponentialBinProbabilities(double sigma_total, double L,
                                     const std::vector<double>& edges);

enum class Target : std::uint8_t { Hydrogen, Carbon };

// One nuclide of the medium with its cross sections at the incident energy.
struct MediumComponent {
    Target target = Target::Hydrogen;
    double number_density = 0.0; // atoms per mm^3
    double sigma_elastic = 0.0;  // MT2 [mm^2]
    double sigma_capture = 0.0;  // MT102 [mm^2]
    double sigma_total = 0.0;    // MT1 [mm^2]
};

inline constexpr std::size_t kFirstInteractionCategories = 5;

// Probabilities of the first hadronic interaction of the primary neutron,
// indexed by the schema::FirstInteraction code (Unknown = no interaction).
// The categories follow the classification of the simulation: elastic on H
// and on C, non-elastic on C without capture (sigma_total - sigma_elastic -
// sigma_capture), and Other for captures and any non-elastic channel on H.
struct FirstInteractionProbabilities {
    double sigma_total = 0.0; // macroscopic, sum_i n_i sigma_total,i [1/mm]
    // P(category); entry Unknown is P(no interaction) = e^{-Sigma L}.
    std::array<double, kFirstInteractionCategories> unconditional{};
    // P(category | an interaction); entry Unknown is 0.
    std::array<double, kFirstInteractionCategories> conditional{};

    [[nodiscard]] double probability(schema::FirstInteraction category) const {
        return unconditional.at(static_cast<std::size_t>(category));
    }
    [[nodiscard]] double probabilityGivenInteraction(schema::FirstInteraction category) const {
        return conditional.at(static_cast<std::size_t>(category));
    }
};

// First-interaction probabilities for a path L [mm] through the medium.
// Throws std::invalid_argument for an empty medium, negative densities or
// cross sections, sigma_elastic + sigma_capture above sigma_total beyond
// rounding, a zero total macroscopic cross section, or L <= 0.
[[nodiscard]] FirstInteractionProbabilities
firstInteraction(const std::vector<MediumComponent>& components, double L);

} // namespace riptide::sim::checks
