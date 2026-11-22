#pragma once

// From the simulated truth to what the detector would report: the HP
// target-mass convention is removed from the proton energy, then the energy
// and angular resolutions are applied.
//
// Two measurement models, as in the toy model:
//  - case A uses the recoil angle theta_p alone, smeared with N(0, sigma_theta)
//    in one dimension, the model its likelihood assumes;
//  - cases B and C use the full track direction, displaced by an isotropic
//    Gaussian in the tangent plane with sigma_theta per component, so that its
//    projection on any meridian plane is again N(0, sigma_theta).
// The two use independent draws; deriving theta_p from the smeared direction
// is a later refinement.
//
// Every random function takes the engine as an argument. applyResolution
// gives each event its own engine, seeded from (seed, event_id): the noise of
// an event does not depend on which other events were selected or in what
// order they were read.

#include <cstdint>
#include <utility>

#include "riptide/common/event_sample.hpp"
#include "riptide/common/random.hpp"

namespace riptide::io {

struct Direction {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Resolution {
    double sigma_Ep = 0.0;    // MeV
    double sigma_theta = 0.0; // rad, per component for directions
};

// Throws std::invalid_argument unless both widths are finite and >= 0.
void requireValidResolution(const Resolution& resolution);

// Two unit vectors (e1, e2) such that (e1, e2, axis) is a right-handed
// orthonormal frame; `axis` must be a unit vector.
[[nodiscard]] std::pair<Direction, Direction> perpendicularBasis(const Direction& axis) noexcept;

// Unit vector at polar angle theta and azimuth phi about the unit vector `axis`.
[[nodiscard]] Direction directionAroundAxis(const Direction& axis, double theta,
                                            double phi) noexcept;

// Ep + N(0, sigma_Ep).
[[nodiscard]] double smearProtonEnergy(double Ep, double sigma_Ep, random::Engine& engine);

// theta_p + N(0, sigma_theta). The result is not folded into [0, pi]: the
// case-A likelihood models exactly this Gaussian.
[[nodiscard]] double smearRecoilAngle(double theta_p, double sigma_theta, random::Engine& engine);

// Whether an event has a proton track. The event files write a zero direction
// when the first n-p produced no proton (truth_valid_np = 0), so a zero vector
// means "no track" and a unit vector a track. Anything else is a damaged input
// and throws std::invalid_argument: smearing it would give NaN silently.
[[nodiscard]] bool hasTrack(const Direction& track);

// The track moved by a Gaussian offset in the tangent plane, two draws of
// N(0, sigma_theta); `track` must be a unit vector, the result is one.
[[nodiscard]] Direction smearDirection(const Direction& track, double sigma_theta,
                                       random::Engine& engine);

// Fills truth_Ep_phys from truth_Ep and truth_En with the HP target mass of
// the run. Throws std::invalid_argument unless the mass is positive and finite
// and the two input columns have one entry per event.
void convertHpEnergies(EventSample& sample, double target_mass);

// Fills reco_Ep (from truth_Ep_phys), reco_theta_p and reco_p_dir_* with the
// per-event engine of (seed, event_id). Draw order within an event: Ep, theta_p,
// then the two direction offsets. An event without a track (see hasTrack) gets
// reco_p_dir = (0, 0, 0) and no direction draws; its Ep and theta_p draws are
// unchanged. Throws std::invalid_argument on an invalid resolution, a negative
// event_id, a truth direction that is neither unit nor zero, or missing input
// columns.
void applyResolution(EventSample& sample, const Resolution& resolution, std::uint64_t seed);

} // namespace riptide::io
