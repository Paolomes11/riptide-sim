#include "riptide/io/response.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "riptide/common/hp_kinematics.hpp"

namespace riptide::io {

namespace {

// Above this |axis.x| the helper axis switches from x to y, so that the cross
// product with the helper never comes close to zero.
constexpr double kHelperSwitch = 0.9;

Direction cross(const Direction& a, const Direction& b) noexcept {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

Direction normalized(const Direction& v) noexcept {
    const double norm = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return {v.x / norm, v.y / norm, v.z / norm};
}

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

std::pair<Direction, Direction> perpendicularBasis(const Direction& axis) noexcept {
    const Direction helper =
        std::abs(axis.x) < kHelperSwitch ? Direction{1.0, 0.0, 0.0} : Direction{0.0, 1.0, 0.0};
    const Direction e1 = normalized(cross(axis, helper));
    return {e1, cross(axis, e1)};
}

Direction directionAroundAxis(const Direction& axis, double theta, double phi) noexcept {
    const auto [e1, e2] = perpendicularBasis(axis);
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);
    const double along_e1 = sin_theta * std::cos(phi);
    const double along_e2 = sin_theta * std::sin(phi);
    return {cos_theta * axis.x + along_e1 * e1.x + along_e2 * e2.x,
            cos_theta * axis.y + along_e1 * e1.y + along_e2 * e2.y,
            cos_theta * axis.z + along_e1 * e1.z + along_e2 * e2.z};
}

double smearProtonEnergy(double Ep, double sigma_Ep, random::Engine& engine) {
    return Ep + sigma_Ep * random::standardNormal(engine);
}

double smearRecoilAngle(double theta_p, double sigma_theta, random::Engine& engine) {
    return theta_p + sigma_theta * random::standardNormal(engine);
}

Direction smearDirection(const Direction& track, double sigma_theta, random::Engine& engine) {
    // Two draws in a fixed order (first e1, then e2): with a function call per
    // argument the order of evaluation would be unspecified.
    const double offset_1 = sigma_theta * random::standardNormal(engine);
    const double offset_2 = sigma_theta * random::standardNormal(engine);
    return directionAroundAxis(track, std::hypot(offset_1, offset_2),
                               std::atan2(offset_2, offset_1));
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
        const auto reco = smearDirection(
            {sample.truth_p_dir_x[i], sample.truth_p_dir_y[i], sample.truth_p_dir_z[i]},
            resolution.sigma_theta, engine);
        sample.reco_p_dir_x[i] = reco.x;
        sample.reco_p_dir_y[i] = reco.y;
        sample.reco_p_dir_z[i] = reco.z;
    }
}

} // namespace riptide::io
