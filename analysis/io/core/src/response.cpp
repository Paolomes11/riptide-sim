#include "riptide/io/response.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "riptide/common/hp_kinematics.hpp"

namespace riptide::io {

namespace {

// Tolerance on the norm of a truth direction: the simulation writes unit
// vectors to about 1e-16 and exact zeros, so 1e-9 separates the two cases with
// a wide margin while still catching a damaged column.
constexpr double kTrackNormTolerance = 1e-9;

void requireColumn(const std::vector<double>& column, std::size_t n, const char* name) {
    if (column.size() != n) {
        throw std::invalid_argument(std::string{"column "} + name + " has not one entry per event");
    }
}

} // namespace

void requireValidResolution(const Resolution& resolution) {
    const auto valid = [](double sigma) { return std::isfinite(sigma) && sigma >= 0.0; };
    if (!valid(resolution.sigma_Ep) || !valid(resolution.sigma_theta)) {
        throw std::invalid_argument("resolutions must be finite and >= 0");
    }
}

double smearProtonEnergy(double Ep, double sigma_Ep, random::Engine& engine) {
    return Ep + sigma_Ep * random::standardNormal(engine);
}

double smearRecoilAngle(double theta_p, double sigma_theta, random::Engine& engine) {
    return theta_p + sigma_theta * random::standardNormal(engine);
}

bool hasTrack(const Direction& track) {
    const double norm = std::sqrt(track.x * track.x + track.y * track.y + track.z * track.z);
    if (std::abs(norm - 1.0) < kTrackNormTolerance) {
        return true;
    }
    if (norm < kTrackNormTolerance) {
        return false;
    }
    // Also reached by NaN, for which both comparisons above are false.
    throw std::invalid_argument("truth proton direction is neither a unit vector nor zero");
}

void convertHpEnergies(EventSample& sample, double target_mass) {
    if (!std::isfinite(target_mass) || target_mass <= 0.0) {
        throw std::invalid_argument("HP target mass must be positive and finite");
    }
    const auto n = sample.size();
    requireColumn(sample.truth_Ep, n, "truth_Ep");
    requireColumn(sample.truth_En, n, "truth_En");
    sample.truth_Ep_phys.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        sample.truth_Ep_phys[i] = kinematics::physicalProtonEnergyFromHp(
            sample.truth_Ep[i], sample.truth_En[i], target_mass);
    }
}

void applyResolution(EventSample& sample, const Resolution& resolution, std::uint64_t seed) {
    requireValidResolution(resolution);
    const auto n = sample.size();
    requireColumn(sample.truth_Ep_phys, n, "truth_Ep_phys");
    requireColumn(sample.truth_theta_p, n, "truth_theta_p");
    requireColumn(sample.truth_p_dir_x, n, "truth_p_dir_x");
    requireColumn(sample.truth_p_dir_y, n, "truth_p_dir_y");
    requireColumn(sample.truth_p_dir_z, n, "truth_p_dir_z");
    sample.reco_Ep.resize(n);
    sample.reco_theta_p.resize(n);
    sample.reco_p_dir_x.resize(n);
    sample.reco_p_dir_y.resize(n);
    sample.reco_p_dir_z.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (sample.event_id[i] < 0) {
            throw std::invalid_argument("negative event_id");
        }
        auto engine = random::perEventEngine(seed, static_cast<std::uint64_t>(sample.event_id[i]));
        sample.reco_Ep[i] = smearProtonEnergy(sample.truth_Ep_phys[i], resolution.sigma_Ep, engine);
        sample.reco_theta_p[i] =
            smearRecoilAngle(sample.truth_theta_p[i], resolution.sigma_theta, engine);
        const Direction track{sample.truth_p_dir_x[i], sample.truth_p_dir_y[i],
                              sample.truth_p_dir_z[i]};
        // No track, no direction to smear: the zero vector stays the marker.
        const Direction reco =
            hasTrack(track) ? smearDirection(track, resolution.sigma_theta, engine) : Direction{};
        sample.reco_p_dir_x[i] = reco.x;
        sample.reco_p_dir_y[i] = reco.y;
        sample.reco_p_dir_z[i] = reco.z;
    }
}

} // namespace riptide::io
