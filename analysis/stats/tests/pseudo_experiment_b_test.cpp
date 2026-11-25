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
#include "riptide/stats/acceptance.hpp"
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

namespace {

constexpr double kHalfPi = 0.5 * pi;
constexpr std::size_t kAccThetaBins = 18;
constexpr std::size_t kAccEvents = 200000;

// A ratio with a forward deficit (half the isotropic rate at 0 deg) that
// deepens with E_n, on 1 MeV bins from 10 to 20 MeV.
double accShape(std::size_t k, std::size_t t) {
    const double depth = 0.4 + (0.02 * static_cast<double>(k));
    return t < 6 ? 1.0 - (depth * (1.0 - (static_cast<double>(t) / 6.0))) : 1.0;
}

double isoShape(std::size_t /*k*/, std::size_t /*t*/) {
    return 1.0;
}

st::AcceptanceRatio accRatio(double (*shape)(std::size_t, std::size_t) = accShape) {
    st::AcceptanceHistogram h;
    for (std::size_t k = 0; k <= 10; ++k) {
        h.En_edges.push_back(10.0 + static_cast<double>(k));
    }
    for (std::size_t t = 0; t <= kAccThetaBins; ++t) {
        h.theta_edges.push_back(kHalfPi * static_cast<double>(t) /
                                static_cast<double>(kAccThetaBins));
    }
    h.theta_edges.back() = kHalfPi;
    for (std::size_t k = 0; k < 10; ++k) {
        const double En = 10.5 + static_cast<double>(k);
        for (std::size_t t = 0; t < kAccThetaBins; ++t) {
            h.counts.push_back(
                shape(k, t) * 1e4 *
                st::isotropicRecoilAngleMass(En, h.theta_edges[t], h.theta_edges[t + 1]));
        }
    }
    return st::AcceptanceRatio(h);
}

std::size_t thetaBin(double theta) {
    return std::min(kAccThetaBins - 1, static_cast<std::size_t>(theta / kHalfPi * kAccThetaBins));
}

// Expected theta_p bin probabilities of q(theta | En) = r pi_iso / norm, by
// a fine midpoint rule in u = Ep / Ep_max (pi_iso is uniform in u).
std::vector<double> expectedThetaBins(const st::AcceptanceRatio& ratio, double En) {
    constexpr std::size_t kSteps = 200000;
    std::vector<double> p(kAccThetaBins, 0.0);
    double total = 0.0;
    for (std::size_t s = 0; s < kSteps; ++s) {
        const double u = (static_cast<double>(s) + 0.5) / static_cast<double>(kSteps);
        const double theta = kin::protonAngleFromEnergy(u * kin::maxProtonEnergy(En), En);
        const double w = ratio(En, theta);
        p[thetaBin(theta)] += w;
        total += w;
    }
    for (double& value : p) {
        value /= total;
    }
    return p;
}

double thetaChiSquare(const std::vector<double>& theta, const std::vector<double>& p) {
    std::vector<double> counts(kAccThetaBins, 0.0);
    for (const double value : theta) {
        counts[thetaBin(value)] += 1.0;
    }
    const auto n = static_cast<double>(theta.size());
    double chi2 = 0.0;
    for (std::size_t t = 0; t < kAccThetaBins; ++t) {
        const double expected = n * p[t];
        chi2 += (counts[t] - expected) * (counts[t] - expected) / expected;
    }
    return chi2;
}

// chi^2 with 17 degrees of freedom: P(chi^2 > 45) ~ 2e-4. The seed is fixed,
// so the check is deterministic; the threshold only says how unlucky a
// correct sampler would have to be.
constexpr double kChi2Limit17 = 45.0;

} // namespace

TEST_CASE("the accepted Case B toy follows q(theta | E_n)", "[stats][pseudo_experiment_b]") {
    const st::AcceptanceRatio ratio = accRatio();
    for (const double En : {11.0, 14.0, 19.0}) {
        // A 1 keV wide domain fixes E_n: q changes far less than the
        // statistical error of the bins over that width.
        const st::EnergyDomain narrow{En, En + 1e-3, 2};
        const auto sample =
            st::caseBPseudoSample(narrow, kOmega, kAccEvents, kToyResolution, kSeed, ratio);
        REQUIRE(sample.size() == kAccEvents);
        const double En_mid = En + 0.5e-3;
        const double chi2 = thetaChiSquare(sample.truth_theta_p, expectedThetaBins(ratio, En_mid));
        INFO("E_n = " << En << " chi2 = " << chi2);
        CHECK(chi2 < kChi2Limit17);
        // The isotropic distribution is clearly rejected: the test has power.
        const double chi2_iso =
            thetaChiSquare(sample.truth_theta_p, expectedThetaBins(accRatio(isoShape), En_mid));
        CHECK(chi2_iso > 10.0 * kChi2Limit17);
        // The kept event is a full Case B event about omega.
        for (std::size_t k = 0; k < sample.size(); ++k) {
            REQUIRE_THAT(geo::angleBetween(sample.truth_track[k], kOmega),
                         WithinAbs(sample.truth_theta_p[k], 1e-9));
        }
    }
}

TEST_CASE("the accepted Case B toy keeps E_n and the azimuth flat and is reproducible",
          "[stats][pseudo_experiment_b]") {
    const st::AcceptanceRatio ratio = accRatio();
    constexpr std::size_t kN = 20000;
    const auto plain = st::caseBPseudoSample(kDomain, kOmega, kN, kToyResolution, kSeed);
    const auto a = st::caseBPseudoSample(kDomain, kOmega, kN, kToyResolution, kSeed, ratio);
    const auto b = st::caseBPseudoSample(kDomain, kOmega, kN, kToyResolution, kSeed, ratio);
    for (std::size_t k = 0; k < kN; ++k) {
        // E_n is the first draw of each engine, before any accept-reject.
        REQUIRE_THAT(a.truth_En[k], WithinULP(plain.truth_En[k], 0));
    }
    requireSameEvents(a, 0, b, kN);
    constexpr std::size_t kTail = 15000;
    const auto tail =
        st::caseBPseudoSample(kDomain, kOmega, kN - kTail, kToyResolution, kSeed, ratio, kTail);
    requireSameEvents(a, kTail, tail, tail.size());

    const auto [e1, e2] = geo::perpendicularBasis(kOmega);
    std::vector<double> u(kN);
    for (std::size_t k = 0; k < kN; ++k) {
        const auto& t = a.truth_track[k];
        u[k] = (std::atan2(geo::dot(t, e2), geo::dot(t, e1)) + pi) / (2.0 * pi);
    }
    CHECK(ksDistance(u, uniformCdf) < ksCritical(kN));

    // r = 1 keeps the first try: the truth is the plain one. The noise
    // comes after the accept draw, so only the truth is compared.
    const auto flat =
        st::caseBPseudoSample(kDomain, kOmega, kN, kToyResolution, kSeed, accRatio(isoShape));
    for (std::size_t k = 0; k < kN; ++k) {
        REQUIRE_THAT(flat.truth_theta_p[k], WithinULP(plain.truth_theta_p[k], 0));
        REQUIRE_THAT(flat.truth_track[k].x, WithinULP(plain.truth_track[k].x, 0));
    }
    // ... and its distributions are the ideal ones.
    std::vector<double> angle(kN);
    std::vector<double> delta(kN);
    for (std::size_t k = 0; k < kN; ++k) {
        angle[k] = geo::angleBetween(flat.truth_track[k], kOmega);
        delta[k] = geo::angleBetween(flat.reco_track[k], flat.truth_track[k]);
    }
    CHECK(ksDistance(angle, thetaCdf) < ksCritical(kN));
    const double s = kToyResolution.sigma_theta;
    const auto rayleigh_cdf = [s](double d) { return 1.0 - std::exp(-d * d / (2.0 * s * s)); };
    CHECK(ksDistance(delta, rayleigh_cdf) < ksCritical(kN));

    CHECK_THROWS_AS(st::caseBPseudoSample(kDomain, kOmega, 1, {0.0, 0.08}, kSeed, ratio),
                    std::invalid_argument);
}
