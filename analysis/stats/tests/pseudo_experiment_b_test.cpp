#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <vector>

#include "riptide/common/direction.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"
#include "riptide/stats/calibration.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"

namespace st = riptide::stats;
namespace kin = riptide::kinematics;
namespace geo = riptide::geometry;
namespace rnd = riptide::random;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;
using riptide::units::pi;

namespace {

constexpr std::uint64_t kSeed = 20260907;
constexpr st::CaseBResolution kToyResolution{0.10, 0.08};
const st::EnergyDomain kDomain{};
const geo::Direction kOmega = geo::normalized({3.0, 2.0, 1.0});

// Every distribution test below draws this many events and accepts the
// one-sample KS distance below the asymptotic critical value at p = 0.01,
// 1.628 / sqrt(n) (about 0.0115). The seed is fixed, so each check is
// deterministic; a failure means a wrong distribution, not bad luck.
constexpr std::size_t kEvents = 20000;
constexpr double kKsCritical001 = 1.628;
// The model CDF of theta_p marginalised over E_n by trapezoids on this many
// points: the integrand is smooth in E_n and nearly constant, so the
// quadrature error is far below the KS critical value.
constexpr std::size_t kQuadraturePoints = 201;

double ksDistance(std::vector<double> values, const std::function<double(double)>& cdf) {
    std::sort(values.begin(), values.end());
    const auto n = static_cast<double>(values.size());
    double distance = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double f = cdf(values[i]);
        const auto index = static_cast<double>(i);
        distance = std::max({distance, std::abs(f - index / n), std::abs(f - (index + 1.0) / n)});
    }
    return distance;
}

double ksCritical(std::size_t n) {
    return kKsCritical001 / std::sqrt(static_cast<double>(n));
}

// P(theta_p <= t) for E_n uniform on the domain: at fixed E_n the recoil Ep
// is uniform on [0, Ep_max], and Ep decreases with theta_p, so
// P(theta_p <= t | E_n) = 1 - Ep(E_n, t) / Ep_max(E_n).
double thetaCdf(double t) {
    const double step =
        (kDomain.En_max - kDomain.En_min) / static_cast<double>(kQuadraturePoints - 1);
    double sum = 0.0;
    for (std::size_t i = 0; i < kQuadraturePoints; ++i) {
        const double En = kDomain.En_min + (static_cast<double>(i) * step);
        const double value = 1.0 - (kin::protonEnergy(En, t) / kin::maxProtonEnergy(En));
        const double weight = (i == 0 || i + 1 == kQuadraturePoints) ? 0.5 : 1.0;
        sum += weight * value;
    }
    return sum / static_cast<double>(kQuadraturePoints - 1);
}

double uniformCdf(double x) {
    return std::clamp(x, 0.0, 1.0);
}

const st::CaseBPseudoSample& toySample() {
    static const st::CaseBPseudoSample sample =
        st::caseBPseudoSample(kDomain, kOmega, kEvents, kToyResolution, kSeed);
    return sample;
}

void requireSameEvents(const st::CaseBPseudoSample& a, std::size_t offset_a,
                       const st::CaseBPseudoSample& b, std::size_t count) {
    for (std::size_t k = 0; k < count; ++k) {
        const std::size_t i = offset_a + k;
        CHECK_THAT(a.truth_En[i], WithinULP(b.truth_En[k], 0));
        CHECK_THAT(a.truth_Ep[i], WithinULP(b.truth_Ep[k], 0));
        CHECK_THAT(a.truth_theta_p[i], WithinULP(b.truth_theta_p[k], 0));
        CHECK_THAT(a.reco_Ep[i], WithinULP(b.reco_Ep[k], 0));
        CHECK_THAT(a.reco_track[i].x, WithinULP(b.reco_track[k].x, 0));
        CHECK_THAT(a.reco_track[i].y, WithinULP(b.reco_track[k].y, 0));
        CHECK_THAT(a.reco_track[i].z, WithinULP(b.reco_track[k].z, 0));
    }
}

} // namespace

TEST_CASE("Case B pseudo-events obey the exact kinematics", "[stats][pseudo_experiment_b]") {
    const auto& sample = toySample();
    REQUIRE(sample.size() == kEvents);
    for (std::size_t k = 0; k < sample.size(); ++k) {
        const double En = sample.truth_En[k];
        CHECK(En >= kDomain.En_min);
        CHECK(En < kDomain.En_max);
        CHECK(sample.truth_Ep[k] >= 0.0);
        CHECK(sample.truth_Ep[k] <= kin::maxProtonEnergy(En));
        CHECK(sample.truth_theta_p[k] >= 0.0);
        CHECK(sample.truth_theta_p[k] <= pi / 2.0);
        // The track lies at theta_p from omega_n; angleBetween is accurate
        // to a few ulp of pi, the composition of cos/sin to a few more.
        CHECK_THAT(geo::angleBetween(sample.truth_track[k], kOmega),
                   WithinAbs(sample.truth_theta_p[k], 1e-12));
        CHECK_THAT(geo::dot(sample.reco_track[k], sample.reco_track[k]), WithinAbs(1.0, 1e-12));
    }
}

TEST_CASE("Case B pseudo-events follow the model distributions", "[stats][pseudo_experiment_b]") {
    const auto& sample = toySample();
    const double critical = ksCritical(sample.size());

    SECTION("E_n is uniform on the domain") {
        std::vector<double> u(sample.size());
        for (std::size_t k = 0; k < u.size(); ++k) {
            u[k] = (sample.truth_En[k] - kDomain.En_min) / (kDomain.En_max - kDomain.En_min);
        }
        CHECK(ksDistance(u, uniformCdf) < critical);
    }
    SECTION("Ep / Ep_max(E_n) is uniform on [0, 1]") {
        std::vector<double> u(sample.size());
        for (std::size_t k = 0; k < u.size(); ++k) {
            u[k] = sample.truth_Ep[k] / kin::maxProtonEnergy(sample.truth_En[k]);
        }
        CHECK(ksDistance(u, uniformCdf) < critical);
    }
    SECTION("the track-omega angle follows the E_n-marginal CDF") {
        std::vector<double> angle(sample.size());
        for (std::size_t k = 0; k < angle.size(); ++k) {
            angle[k] = geo::angleBetween(sample.truth_track[k], kOmega);
        }
        CHECK(ksDistance(angle, thetaCdf) < critical);
    }
    SECTION("the azimuth about omega is uniform") {
        const auto [e1, e2] = geo::perpendicularBasis(kOmega);
        std::vector<double> u(sample.size());
        for (std::size_t k = 0; k < u.size(); ++k) {
            const auto& t = sample.truth_track[k];
            u[k] = (std::atan2(geo::dot(t, e2), geo::dot(t, e1)) + pi) / (2.0 * pi);
        }
        CHECK(ksDistance(u, uniformCdf) < critical);
    }
    SECTION("the Ep noise is N(0, sigma_Ep)") {
        std::vector<double> z(sample.size());
        for (std::size_t k = 0; k < z.size(); ++k) {
            z[k] = (sample.reco_Ep[k] - sample.truth_Ep[k]) / kToyResolution.sigma_Ep;
        }
        const auto normal_cdf = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };
        CHECK(ksDistance(z, normal_cdf) < critical);
    }
    SECTION("the track offset angle is Rayleigh(sigma_theta)") {
        // Two independent N(0, sigma) offsets in the tangent plane: their
        // norm, which is the angle moved, has CDF 1 - exp(-d^2 / 2 sigma^2).
        std::vector<double> delta(sample.size());
        for (std::size_t k = 0; k < delta.size(); ++k) {
            delta[k] = geo::angleBetween(sample.reco_track[k], sample.truth_track[k]);
        }
        const double s = kToyResolution.sigma_theta;
        const auto rayleigh_cdf = [s](double d) { return 1.0 - std::exp(-d * d / (2.0 * s * s)); };
        CHECK(ksDistance(delta, rayleigh_cdf) < critical);
    }
}

TEST_CASE("isotropicDirection is uniform on the sphere", "[stats][pseudo_experiment_b]") {
    std::vector<double> u_cos(kEvents);
    std::vector<double> u_phi(kEvents);
    for (std::size_t i = 0; i < kEvents; ++i) {
        auto engine = rnd::perEventEngine(kSeed, i);
        const geo::Direction d = st::isotropicDirection(engine);
        CHECK_THAT(geo::dot(d, d), WithinAbs(1.0, 1e-12));
        u_cos[i] = (d.z + 1.0) / 2.0;
        u_phi[i] = (std::atan2(d.y, d.x) + pi) / (2.0 * pi);
    }
    const double critical = ksCritical(kEvents);
    CHECK(ksDistance(u_cos, uniformCdf) < critical);
    CHECK(ksDistance(u_phi, uniformCdf) < critical);
}

TEST_CASE("Case B pseudo-events do not depend on the sample size", "[stats][pseudo_experiment_b]") {
    constexpr std::uint64_t kFirst = 1000;
    const auto small = st::caseBPseudoSample(kDomain, kOmega, 5, kToyResolution, kSeed, kFirst);
    const auto large = st::caseBPseudoSample(kDomain, kOmega, 10, kToyResolution, kSeed, kFirst);
    REQUIRE(small.size() == 5);
    REQUIRE(large.size() == 10);
    requireSameEvents(large, 0, small, small.size());
    // Event k uses id first + k: a sample starting two ids later is the tail.
    const auto shifted =
        st::caseBPseudoSample(kDomain, kOmega, 8, kToyResolution, kSeed, kFirst + 2);
    requireSameEvents(large, 2, shifted, shifted.size());
}

TEST_CASE("Case B pseudo-experiments use the reserved truth id", "[stats][pseudo_experiment_b]") {
    constexpr std::uint64_t kSection = 3;
    constexpr std::uint64_t kIndex = 7;
    constexpr std::size_t kN = 20;
    const auto experiment =
        st::caseBPseudoExperiment(kDomain, kN, kToyResolution, kSeed, kSection, kIndex);
    const std::uint64_t first = st::experimentFirstEventId(kSection, kIndex);
    auto engine = rnd::perEventEngine(kSeed, first + st::kMaxExperimentEvents);
    const geo::Direction truth = st::isotropicDirection(engine);
    CHECK_THAT(experiment.truth_omega_n.x, WithinULP(truth.x, 0));
    CHECK_THAT(experiment.truth_omega_n.y, WithinULP(truth.y, 0));
    CHECK_THAT(experiment.truth_omega_n.z, WithinULP(truth.z, 0));
    REQUIRE(experiment.sample.size() == kN);
    const auto expected = st::caseBPseudoSample(kDomain, truth, kN, kToyResolution, kSeed, first);
    requireSameEvents(experiment.sample, 0, expected, kN);

    const auto at_direction = st::caseBPseudoExperimentAtDirection(
        truth, kDomain, kN, kToyResolution, kSeed, kSection, kIndex);
    requireSameEvents(at_direction.sample, 0, expected, kN);

    // Another experiment draws another truth.
    const auto other =
        st::caseBPseudoExperiment(kDomain, kN, kToyResolution, kSeed, kSection, kIndex + 1);
    CHECK(geo::angleBetween(other.truth_omega_n, experiment.truth_omega_n) > 1e-6);
}

TEST_CASE("Case B pseudo-experiments reject invalid input", "[stats][pseudo_experiment_b]") {
    CHECK_THROWS_AS(st::caseBPseudoSample(kDomain, kOmega, 1, {0.0, 0.08}, kSeed),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::caseBPseudoSample(kDomain, kOmega, 1, {0.10, -0.08}, kSeed),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::caseBPseudoSample({0.0, 20.0, 10}, kOmega, 1, kToyResolution, kSeed),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::caseBPseudoSample({20.0, 10.0, 10}, kOmega, 1, kToyResolution, kSeed),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::caseBPseudoExperiment(kDomain, st::kMaxExperimentEvents + 1, kToyResolution,
                                              kSeed, 0, 0),
                    std::invalid_argument);
    CHECK_THROWS_AS(st::caseBPseudoExperimentAtDirection(
                        kOmega, kDomain, st::kMaxExperimentEvents + 1, kToyResolution, kSeed, 0, 0),
                    std::invalid_argument);
    constexpr std::uint64_t kSectionOutOfRange = std::uint64_t{1} << 40;
    CHECK_THROWS_AS(
        st::caseBPseudoExperiment(kDomain, 1, kToyResolution, kSeed, kSectionOutOfRange, 0),
        std::invalid_argument);
    CHECK_NOTHROW(st::caseBPseudoSample(kDomain, kOmega, 0, kToyResolution, kSeed));
}
