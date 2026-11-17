#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

namespace kin = riptide::kinematics;
namespace u = riptide::units;

struct LabRecoil {
    double kinetic = 0.0;
    double theta = 0.0;
};

// Independent route to the recoil: build the CM state from the invariant mass,
// scatter the projectile by cos_cm, boost the target back to the lab. Shares
// no formula with the header except the masses.
LabRecoil recoilByBoost(double t_projectile, double cos_cm, double m_projectile, double m_target) {
    const double e_projectile = t_projectile + m_projectile;
    const double p_lab = std::sqrt(t_projectile * (t_projectile + 2.0 * m_projectile));
    const double s_invariant =
        m_projectile * m_projectile + m_target * m_target + 2.0 * m_target * e_projectile;
    const double sqrt_s = std::sqrt(s_invariant);
    const double beta = p_lab / (e_projectile + m_target);
    const double gamma = (e_projectile + m_target) / sqrt_s;
    const double p_cm = m_target * p_lab / sqrt_s;
    const double e_target_cm = std::sqrt(p_cm * p_cm + m_target * m_target);
    // The target recoils opposite to the scattered projectile in the CM.
    const double pz_cm = -p_cm * cos_cm;
    const double pt = p_cm * std::sqrt(1.0 - cos_cm * cos_cm);
    const double pz_lab = gamma * (pz_cm + beta * e_target_cm);
    const double e_lab = gamma * (e_target_cm + beta * pz_cm);
    return {e_lab - m_target, std::atan2(pt, pz_lab)};
}

} // namespace

TEST_CASE("toy oracle: En = 20 MeV, theta_p = 45 deg gives Ep = 9.947 MeV", "[kinematics]") {
    // (a) value reported by the toy model with equal masses m = m_p.
    constexpr double kToyEp = 9.947 * u::MeV;
    constexpr double kToyRounding = 0.5e-3 * u::MeV;
    const double theta_p = 45.0 * u::deg;
    CHECK_THAT(kin::recoilEnergy(20.0 * u::MeV, theta_p, kin::kProtonMass, kin::kProtonMass),
               WithinAbs(kToyEp, kToyRounding));
    // m_n != m_p moves Ep by ~5e-6 MeV only, still within the toy rounding.
    CHECK_THAT(kin::protonEnergy(20.0 * u::MeV, theta_p), WithinAbs(kToyEp, kToyRounding));
}

TEST_CASE("equal masses reduce to the toy closed form", "[kinematics]") {
    const double En = GENERATE(10.0, 15.0, 20.0);
    const double theta_p = GENERATE(0.0, 0.3, 0.785, 1.2, 1.5);
    const double m = kin::kProtonMass;
    const double sin_theta = std::sin(theta_p);
    const double cos_theta = std::cos(theta_p);
    const double toy =
        2.0 * m * En * cos_theta * cos_theta / (2.0 * m + En * sin_theta * sin_theta);
    CHECK_THAT(kin::recoilEnergy(En, theta_p, m, m), WithinRel(toy, 1e-13));
}

TEST_CASE("forward recoil carries the maximum energy", "[kinematics]") {
    const double En = GENERATE(10.0, 14.0, 20.0);
    CHECK_THAT(kin::protonEnergy(En, 0.0), WithinRel(kin::maxProtonEnergy(En), 1e-14));
    // Equal masses: the whole kinetic energy is transferred head-on.
    CHECK_THAT(kin::maxRecoilEnergy(En, kin::kProtonMass, kin::kProtonMass), WithinRel(En, 1e-14));
    // m_n != m_p: the transfer falls short of En, by ~5e-7 relative.
    CHECK(kin::maxProtonEnergy(En) < En);
}

TEST_CASE("non-relativistic limit", "[kinematics]") {
    // For t << m the recoil is 4 m1 m2 / (m1 + m2)^2 t cos^2(theta); the
    // relativistic correction is O(t / m) ~ 1e-9 at 1 eV.
    constexpr double kTiny = 1.0e-6 * u::MeV;
    const double m1 = kin::kNeutronMass;
    const double m2 = kin::kProtonMass;
    const double transfer = 4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2));
    const double theta_p = GENERATE(0.0, 0.5, 1.0, 1.4);
    const double cos_theta = std::cos(theta_p);
    CHECK_THAT(kin::protonEnergy(kTiny, theta_p),
               WithinRel(transfer * kTiny * cos_theta * cos_theta, 1e-8));
}

TEST_CASE("CM map agrees with an explicit Lorentz boost", "[kinematics]") {
    const double En = GENERATE(10.0, 19.5, 20.0);
    const double cos_cm = GENERATE(-0.99, -0.5, 0.0, 0.4, 0.95);
    const LabRecoil boosted = recoilByBoost(En, cos_cm, kin::kNeutronMass, kin::kProtonMass);
    // Ep is linear in cos_cm, exactly.
    CHECK_THAT(kin::protonEnergyFromCosCm(En, cos_cm), WithinRel(boosted.kinetic, 1e-12));
    // The lab-angle formula gives the same energy at the boosted angle.
    CHECK_THAT(kin::protonEnergy(En, boosted.theta), WithinRel(boosted.kinetic, 1e-12));
    // And the map inverts.
    CHECK_THAT(kin::cosCmFromProtonEnergy(boosted.kinetic, En), WithinAbs(cos_cm, 1e-12));
}

TEST_CASE("Legendre polynomials: closed forms and end points", "[kinematics][legendre]") {
    const double x = GENERATE(-1.0, -0.7, -0.2, 0.0, 0.33, 0.9, 1.0);
    const double x2 = x * x;
    CHECK_THAT(kin::legendreP(0, x), WithinAbs(1.0, 1e-15));
    CHECK_THAT(kin::legendreP(1, x), WithinAbs(x, 1e-15));
    CHECK_THAT(kin::legendreP(2, x), WithinAbs(0.5 * (3.0 * x2 - 1.0), 1e-14));
    CHECK_THAT(kin::legendreP(3, x), WithinAbs(0.5 * (5.0 * x2 * x - 3.0 * x), 1e-14));
    CHECK_THAT(kin::legendreP(4, x), WithinAbs((35.0 * x2 * x2 - 30.0 * x2 + 3.0) / 8.0, 1e-14));
    CHECK_THAT(kin::legendreP(-1, x), WithinAbs(0.0, 0.0));
    for (int l = 0; l <= 10; ++l) {
        CHECK_THAT(kin::legendreP(l, 1.0), WithinAbs(1.0, 1e-13));
        CHECK_THAT(kin::legendreP(l, -1.0), WithinAbs(l % 2 == 0 ? 1.0 : -1.0, 1e-13));
    }
}

TEST_CASE("Legendre polynomials are orthogonal on -1 <= x <= 1", "[kinematics][legendre]") {
    // Composite Simpson rule: exact for products of degree <= 3, its error
    // grows with the degree and is ~1e-11 for P_6^2 at this step. A wrong
    // recurrence coefficient would show up at the 1e-2 level, so 1e-9 keeps a
    // wide margin on both sides.
    constexpr int kIntervals = 4000;
    constexpr int kMaxDegree = 6;
    const double step = 2.0 / kIntervals;
    for (int l = 0; l <= kMaxDegree; ++l) {
        for (int m = 0; m <= kMaxDegree; ++m) {
            double sum = 0.0;
            for (int i = 0; i <= kIntervals; ++i) {
                const double x = -1.0 + i * step;
                const double weight = (i == 0 || i == kIntervals) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
                sum += weight * kin::legendreP(l, x) * kin::legendreP(m, x);
            }
            const double integral = sum * step / 3.0;
            const double expected = l == m ? 2.0 / (2.0 * l + 1.0) : 0.0;
            INFO("l = " << l << ", m = " << m);
            CHECK_THAT(integral, WithinAbs(expected, 1e-9));
        }
    }
}
