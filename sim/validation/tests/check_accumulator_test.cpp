#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/hp_recoil.hpp"

namespace checks = riptide::sim::checks;
namespace endf = riptide::sim::endf;
namespace kin = riptide::kinematics;
namespace units = riptide::units;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::schema::FirstInteraction;
using riptide::sim::TruthRecord;
using riptide::sim::Vec3;

namespace {

constexpr std::uint64_t kSeed = 20260907; // common seed of the project
constexpr double kEn = 14.0 * units::MeV;
constexpr double kTolerance = 1e-12;

// Densities and cross sections of the order of a plastic scintillator at
// 14 MeV; the tests only need a plausible medium, not the evaluated values.
constexpr double kHydrogenDensity = 5.0e19; // atoms / mm^3
constexpr double kCarbonDensity = 4.0e19;   // atoms / mm^3
constexpr double kHydrogenElastic = 0.7 * units::barn;
constexpr double kCarbonElastic = 0.8 * units::barn;
constexpr double kCarbonCapture = 1.0e-5 * units::barn; // rare on purpose: it gets folded
constexpr double kCarbonTotal = 1.3 * units::barn;

// Energy-independent table over a range that covers every test energy.
endf::Tab1 constantTable(double value) {
    return endf::Tab1{{2}, {endf::InterpolationLaw::LinLin}, {1.0, 30.0}, {value, value}};
}

checks::ReferenceModel model(const std::string& library, const std::vector<double>& a) {
    checks::ReferenceModel m;
    m.library = library;
    m.np_angular = endf::LegendreTable{
        endf::Frame::CenterOfMass, {2}, {endf::InterpolationLaw::LinLin}, {1.0, 30.0}, {a, a}};
    m.nuclides.push_back(checks::Nuclide{checks::Target::Hydrogen, "H-1", kHydrogenDensity,
                                         constantTable(kHydrogenElastic),
                                         constantTable(kHydrogenElastic), constantTable(0.0)});
    m.nuclides.push_back(checks::Nuclide{checks::Target::Carbon, "C-0", kCarbonDensity,
                                         constantTable(kCarbonTotal), constantTable(kCarbonElastic),
                                         constantTable(kCarbonCapture)});
    return m;
}

const std::vector<double> kCoefficients{0.1, 0.05};

// Pencil beam along +z from outside a 100 mm cube: L = 100 mm.
riptide::sim::RunConfig cubeConfig() {
    riptide::sim::RunConfig config;
    config.geometry.scint_size = {100.0, 100.0, 100.0};
    config.source.En_min = kEn;
    config.source.En_max = kEn;
    config.source.direction = {0.0, 0.0, 1.0};
    config.source.position = {0.0, 0.0, -100.0};
    return config;
}

// With the bare proton as the HP target the HP kinematics is the physical
// one, so the synthetic events below can be built from np_kinematics.
checks::CheckSetup cubeSetup(double target_mass = kin::kProtonMass) {
    return checks::makeCheckSetup(cubeConfig(), target_mass, checks::kDefaultBins,
                                  checks::kDefaultBins);
}

// Target mass of H-1 in the G4NDL 4.7.1 elastic file, in neutron masses.
constexpr double kAwrH1 = 0.999167;

// Lab angle of a recoil proton of energy Ep, inverting recoilEnergy:
// cos^2(theta) = Ep (E + m_p)^2 / (p^2 (2 m_p + Ep)), E the neutron total energy.
double protonAngle(double En, double Ep) {
    const double p_squared = kin::momentumSquared(En, kin::kNeutronMass);
    const double total_plus_target = En + kin::kNeutronMass + kin::kProtonMass;
    const double cos2 =
        Ep * total_plus_target * total_plus_target / (p_squared * (2.0 * kin::kProtonMass + Ep));
    return std::acos(std::sqrt(std::min(cos2, 1.0)));
}

TruthRecord npRecord(double En, double Ep, double depth) {
    TruthRecord r;
    r.En = En;
    r.omega_n = {0.0, 0.0, 1.0};
    r.r0 = {0.0, 0.0, -50.0 + depth};
    r.first_interaction = FirstInteraction::NpElastic;
    r.n_np_scatters = 1;
    r.Ep = Ep;
    r.theta_p = protonAngle(En, Ep);
    r.valid_np = true;
    return r;
}

// Draws events from the very distributions the accumulator expects: the
// first-interaction category, an exponential depth truncated to L, and for
// n-p a CM cosine from f(mu) by rejection.
class Sampler {
public:
    explicit Sampler(const checks::ReferenceModel& m)
        : probabilities_(checks::firstInteraction(m.mediumAt(kEn), kLength)),
          a_(m.coefficientsAt(kEn)) {}

    TruthRecord next() {
        const double r = uniform_(rng_);
        double cumulative = 0.0;
        std::size_t category = 0;
        for (; category + 1 < checks::kFirstInteractionCategories; ++category) {
            cumulative += probabilities_.unconditional.at(category);
            if (r < cumulative) {
                break;
            }
        }
        const auto kind = static_cast<FirstInteraction>(category);
        if (kind == FirstInteraction::Unknown) {
            TruthRecord none;
            none.En = kEn;
            none.omega_n = {0.0, 0.0, 1.0};
            return none;
        }
        const double sigma = probabilities_.sigma_total;
        const double depth = -std::log1p(-uniform_(rng_) * -std::expm1(-sigma * kLength)) / sigma;
        if (kind != FirstInteraction::NpElastic) {
            TruthRecord other;
            other.En = kEn;
            other.omega_n = {0.0, 0.0, 1.0};
            other.r0 = {0.0, 0.0, -50.0 + depth};
            other.first_interaction = kind;
            return other;
        }
        return npRecord(kEn, kin::protonEnergyFromCosCm(kEn, cosine()), depth);
    }

private:
    static constexpr double kLength = 100.0;

    double density(double mu) const {
        double f = 0.5;
        for (std::size_t l = 1; l <= a_.size(); ++l) {
            f += (2.0 * static_cast<double>(l) + 1.0) / 2.0 * a_.at(l - 1) *
                 kin::legendreP(static_cast<int>(l), mu);
        }
        return f;
    }

    double cosine() {
        double bound = 0.5;
        for (std::size_t l = 1; l <= a_.size(); ++l) {
            bound += (2.0 * static_cast<double>(l) + 1.0) / 2.0 * std::abs(a_.at(l - 1));
        }
        for (;;) {
            const double mu = 2.0 * uniform_(rng_) - 1.0;
            if (uniform_(rng_) * bound < density(mu)) {
                return mu;
            }
        }
    }

    checks::FirstInteractionProbabilities probabilities_;
    std::vector<double> a_;
    std::mt19937_64 rng_{kSeed};
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

} // namespace

TEST_CASE("check setup from the run configuration", "[check_accumulator]") {
    const auto setup = cubeSetup();
    CHECK_THAT(setup.chord.length(), WithinAbs(100.0, kTolerance));
    CHECK_THAT(setup.half_size.z, WithinAbs(50.0, kTolerance));

    auto missing = cubeConfig();
    missing.source.position = {80.0, 0.0, -100.0};
    const double m = kin::kProtonMass;
    CHECK_THROWS_AS(checks::makeCheckSetup(missing, m, 10, 10), std::invalid_argument);
    CHECK_THROWS_AS(checks::makeCheckSetup(cubeConfig(), m, 0, 10), std::invalid_argument);
    CHECK_THROWS_AS(checks::makeCheckSetup(cubeConfig(), 0.0, 10, 10), std::invalid_argument);
    CHECK_THROWS_AS(checks::makeCheckSetup(cubeConfig(), std::nan(""), 10, 10),
                    std::invalid_argument);
}

TEST_CASE("accumulator needs distinct reference models", "[check_accumulator]") {
    CHECK_THROWS_AS(checks::CheckAccumulator(cubeSetup(), {}), std::invalid_argument);
    CHECK_THROWS_AS(
        checks::CheckAccumulator(cubeSetup(), {model("a", kCoefficients), model("a", {})}),
        std::invalid_argument);
}

TEST_CASE("monoenergetic expectations equal the single-energy formulas", "[check_accumulator]") {
    const auto m = model("lib", kCoefficients);
    checks::CheckAccumulator accumulator(cubeSetup(), {m});
    constexpr int kEvents = 1000;
    const double Ep = 0.5 * kin::maxProtonEnergy(kEn);
    for (int i = 0; i < kEvents; ++i) {
        accumulator.add(npRecord(kEn, Ep, 30.0));
    }
    const auto results = accumulator.finish();
    REQUIRE(results.np_events == kEvents);

    const auto& u = results.spectrum_u;
    const auto p_u = checks::protonEnergyBinProbabilities(kCoefficients, u.edges);
    REQUIRE(u.per_library.size() == 1);
    for (std::size_t i = 0; i < p_u.size(); ++i) {
        CHECK_THAT(u.per_library.at(0).expected.at(i), WithinRel(kEvents * p_u.at(i), 1e-10));
    }

    const auto probabilities = checks::firstInteraction(m.mediumAt(kEn), 100.0);
    const auto p_z = checks::truncatedExponentialBinProbabilities(probabilities.sigma_total, 100.0,
                                                                  results.depth_z.edges);
    for (std::size_t i = 0; i < p_z.size(); ++i) {
        CHECK_THAT(results.depth_z.per_library.at(0).expected.at(i),
                   WithinRel(kEvents * p_z.at(i), 1e-10));
    }
    for (std::size_t c = 0; c < checks::kFirstInteractionCategories; ++c) {
        CHECK_THAT(results.fractions.per_library.at(0).expected.at(c),
                   WithinRel(kEvents * probabilities.unconditional.at(c), 1e-10));
    }
}

TEST_CASE("a sample from the expected distributions passes every comparison",
          "[check_accumulator]") {
    const auto m = model("lib", kCoefficients);
    checks::CheckAccumulator accumulator(cubeSetup(), {m, model("flat", {})});
    Sampler sampler(m);
    // Large enough for the isotropic model below to fail clearly.
    constexpr int kEvents = 100000;
    for (int i = 0; i < kEvents; ++i) {
        accumulator.add(sampler.next());
    }
    const auto results = accumulator.finish();
    CHECK(results.events == kEvents);
    CHECK(results.counters.mustBeZeroTotal() == 0);
    CHECK(results.counters.u_overflow == 0);
    CHECK(results.counters.np_first_without_proton == 0);
    REQUIRE(results.np_events > 1000);

    // A p-value below this would be a 1-in-10^4 fluctuation for a correct
    // expectation; the seed is fixed, so the test is deterministic anyway.
    constexpr double kMinP = 1e-4;
    const auto& u = results.spectrum_u.per_library.at(0);
    REQUIRE(u.comparison.has_value());
    CHECK(u.comparison->p_value > kMinP);
    const auto& z = results.depth_z.per_library.at(0);
    REQUIRE(z.comparison.has_value());
    CHECK(z.comparison->p_value > kMinP);
    // Both models share the medium, so the depth agrees for both.
    CHECK(results.depth_z.per_library.at(1).comparison->p_value > kMinP);

    const auto& f = results.fractions.per_library.at(0);
    REQUIRE(f.comparison.has_value());
    CHECK(f.comparison->p_value > kMinP);
    // Other (capture on C) is expected far below kMinExpectedCount: folded.
    CHECK(std::find(f.categories.begin(), f.categories.end(),
                    static_cast<std::int32_t>(FirstInteraction::Other)) == f.categories.end());

    REQUIRE(results.legendre.moments.has_value());
    REQUIRE(results.legendre.moments->values.size() == kCoefficients.size());
    for (const double pull : results.legendre.per_library.at(0).pulls) {
        CHECK(std::abs(pull) < 4.0);
    }
    // The isotropic model is wrong for this sample: its a_1 pull is large.
    CHECK(std::abs(results.legendre.per_library.at(1).pulls.at(0)) > 5.0);
    CHECK(results.spectrum_u.per_library.at(1).comparison->p_value < kMinP);

    CHECK(results.Ep_residuals.samples == results.np_events);
    CHECK(results.Ep_residuals.max_abs < 1e-9 * units::MeV);
}

TEST_CASE("consistency counters", "[check_accumulator]") {
    checks::CheckAccumulator accumulator(cubeSetup(), {model("lib", kCoefficients)});
    const double Ep_max = kin::maxProtonEnergy(kEn);

    auto outside_source = npRecord(15.0, 0.5 * Ep_max, 10.0);
    accumulator.add(outside_source);

    auto tilted = npRecord(kEn, 0.5 * Ep_max, 10.0);
    tilted.omega_n = {0.0, 1e-6, 1.0};
    accumulator.add(tilted);

    auto no_scatter = npRecord(kEn, 0.5 * Ep_max, 10.0);
    no_scatter.n_np_scatters = 0;
    accumulator.add(no_scatter);

    auto off_box = npRecord(kEn, 0.5 * Ep_max, 10.0);
    off_box.r0.x = 51.0;
    accumulator.add(off_box);

    auto zero_energy = npRecord(kEn, 0.5 * Ep_max, 10.0);
    zero_energy.Ep = 0.0;
    accumulator.add(zero_energy);

    auto deuteron = npRecord(kEn, 0.5 * Ep_max, 10.0);
    deuteron.valid_np = false;
    deuteron.Ep = 0.0;
    accumulator.add(deuteron);

    // Ep beyond the kinematic limit: it matches neither kinematics.
    auto overflow = npRecord(kEn, 0.5 * Ep_max, 10.0);
    overflow.Ep = 1.001 * Ep_max;
    accumulator.add(overflow);

    auto backward = npRecord(kEn, 0.5 * Ep_max, 10.0);
    backward.theta_p = 0.5 * units::pi + 1e-3;
    accumulator.add(backward);

    const auto results = accumulator.finish();
    const auto& c = results.counters;
    CHECK(c.En_outside_source == 1);
    CHECK(c.omega_not_beam == 1);
    CHECK(c.valid_np_without_scatter == 1);
    CHECK(c.vertex_outside_box == 1);
    CHECK(c.invalid_proton_energy == 1);
    CHECK(c.theta_p_beyond_half_pi == 1);
    // The overflow and the backward proton fit no kinematics; nothing reaches
    // the u histogram beyond 1, so the overflow counter stays a second guard.
    CHECK(c.np_kinematics_unmatched == 2);
    CHECK(c.u_overflow == 0);
    CHECK(c.mustBeZeroTotal() == 8);
    CHECK(c.np_first_without_proton == 1);
    // With the bare proton as HP target the two kinematics coincide.
    CHECK(c.np_model_ambiguous == 3);

    // The event outside the source range is left out of every comparison.
    CHECK(results.events == 8);
    std::uint64_t compared = 0;
    for (const auto count : results.fractions.observed) {
        compared += count;
    }
    CHECK(compared == 7);
    // Unmatched events are left out of the HP sample.
    double histogrammed = 0.0;
    for (const double o : results.spectrum_u.observed) {
        histogrammed += o;
    }
    CHECK_THAT(histogrammed, WithinAbs(3.0, kTolerance));
    CHECK(results.np_events == 5);
    CHECK(results.np_hp_events == 3);
}

TEST_CASE("events with the HP target mass are modelled to rounding", "[check_accumulator]") {
    const double target_mass = kAwrH1 * kin::kNeutronMass;
    checks::CheckAccumulator hp(cubeSetup(target_mass), {model("lib", kCoefficients)});
    checks::CheckAccumulator bare(cubeSetup(), {model("lib", kCoefficients)});
    // Lab angles from head-on (u = 1, the edge where a wrong mass overflows)
    // to just short of pi/2.
    constexpr int kAngles = 50;
    for (int i = 0; i < kAngles; ++i) {
        const double theta = 0.5 * units::pi * static_cast<double>(i) / kAngles;
        auto r = npRecord(kEn, 0.5 * kin::maxProtonEnergy(kEn), 30.0);
        r.theta_p = theta;
        r.Ep = checks::hpProtonEnergy(kEn, theta, target_mass);
        hp.add(r);
        bare.add(r);
    }
    const auto results = hp.finish();
    CHECK(results.counters.mustBeZeroTotal() == 0);
    CHECK(results.counters.np_model_ambiguous == 0);
    CHECK(results.np_events == kAngles);
    CHECK(results.np_hp_events == kAngles);
    CHECK(results.Ep_residuals.max_abs < 1e-9 * kin::maxProtonEnergy(kEn));
    // The excess over the physical kinematics is the HP convention, ~m_e/m_p.
    const auto& offset = results.Ep_physical_offset;
    CHECK(offset.absolute.samples == kAngles);
    CHECK(offset.relative_min > 5.36e-4);
    CHECK(offset.relative_max < 5.49e-4);
    CHECK(offset.relative_mean > offset.relative_min);
    CHECK(offset.relative_mean < offset.relative_max);
    // Analysed with the bare proton as HP target, the HP events match
    // neither kinematics: a wrong target mass cannot pass unnoticed.
    CHECK(bare.finish().counters.np_kinematics_unmatched == kAngles);
}

TEST_CASE("the rounding-free offset extremes discount a residual at tiny Ep",
          "[check_accumulator]") {
    const double target_mass = kAwrH1 * kin::kNeutronMass;
    checks::CheckAccumulator acc(cubeSetup(target_mass), {model("lib", kCoefficients)});
    // A rounding-sized residual, harmless at MeV but ~1e-10 relative at ~14 eV.
    constexpr double kResidual = 2e-15 * units::MeV;
    constexpr double kForward = 0.3;      // [rad]: Ep of several MeV
    constexpr double kNearGrazing = 1e-3; // [rad] short of pi/2: Ep of ~14 eV
    double exact_min = std::numeric_limits<double>::infinity();
    double exact_max = -std::numeric_limits<double>::infinity();
    for (const double theta : {kForward, 0.5 * units::pi - kNearGrazing}) {
        auto r = npRecord(kEn, 0.5 * kin::maxProtonEnergy(kEn), 30.0);
        r.theta_p = theta;
        const double Ep_hp = checks::hpProtonEnergy(kEn, theta, target_mass);
        const double exact = Ep_hp / kin::protonEnergy(kEn, theta) - 1.0;
        exact_min = std::min(exact_min, exact);
        exact_max = std::max(exact_max, exact);
        r.Ep = Ep_hp + kResidual;
        acc.add(r);
    }
    const auto results = acc.finish();
    REQUIRE(results.np_hp_events == 2);
    CHECK(results.counters.mustBeZeroTotal() == 0);
    const auto& offset = results.Ep_physical_offset;
    // The raw maximum is pushed past the exact one by the residual alone...
    CHECK(offset.relative_max > exact_max + 1e-11);
    // ...while the rounding-free extremes bracket the exact offsets, up to
    // the rounding of the ratio itself (a few ulp of 1).
    constexpr double kRatioRounding = 4.0 * std::numeric_limits<double>::epsilon();
    CHECK(offset.rounding_free_relative_max <= exact_max + kRatioRounding);
    CHECK(offset.rounding_free_relative_min >= exact_min - kRatioRounding);
    CHECK(offset.rounding_free_relative_max <= offset.relative_max);
    CHECK(offset.rounding_free_relative_min >= offset.relative_min);
}

TEST_CASE("non-HP probability is the linear ramp of the model overlap", "[check_accumulator]") {
    CHECK_THAT(checks::nonHpElasticProbability(14.0), WithinAbs(0.0, kTolerance));
    CHECK_THAT(checks::nonHpElasticProbability(checks::kChipsElasticMinEnergy),
               WithinAbs(0.0, kTolerance));
    CHECK_THAT(checks::nonHpElasticProbability(19.75), WithinAbs(0.5, kTolerance));
    CHECK_THAT(checks::nonHpElasticProbability(checks::kHpElasticMaxEnergy),
               WithinAbs(1.0, kTolerance));
    CHECK_THAT(checks::nonHpElasticProbability(25.0), WithinAbs(1.0, kTolerance));
}

namespace {

// Flat source across the HP/CHIPS window with the G4NDL target mass. Each
// event goes to CHIPS (physical kinematics) with probability `ramp(En)`,
// otherwise to HP; the proton angle is uniform below pi/2.
template <typename Ramp> checks::CheckResults windowRun(Ramp ramp, int events) {
    auto config = cubeConfig();
    config.source.En_min = checks::kChipsElasticMinEnergy;
    config.source.En_max = checks::kHpElasticMaxEnergy;
    const double target_mass = kAwrH1 * kin::kNeutronMass;
    checks::CheckAccumulator accumulator(
        checks::makeCheckSetup(config, target_mass, checks::kDefaultBins, checks::kDefaultBins),
        {model("lib", kCoefficients)});
    std::mt19937_64 rng{kSeed};
    std::uniform_real_distribution<double> energy{config.source.En_min, config.source.En_max};
    std::uniform_real_distribution<double> angle{0.0, 0.45 * units::pi};
    std::uniform_real_distribution<double> uniform{0.0, 1.0};
    for (int i = 0; i < events; ++i) {
        const double En = energy(rng);
        const double theta = angle(rng);
        auto r = npRecord(En, 0.5 * kin::maxProtonEnergy(En), 30.0);
        r.theta_p = theta;
        r.Ep = uniform(rng) < ramp(En) ? kin::protonEnergy(En, theta)
                                       : checks::hpProtonEnergy(En, theta, target_mass);
        accumulator.add(r);
    }
    return accumulator.finish();
}

} // namespace

TEST_CASE("events across the HP/CHIPS window are split by their kinematics",
          "[check_accumulator]") {
    constexpr int kEvents = 20000;
    const auto results = windowRun(checks::nonHpElasticProbability, kEvents);
    const auto& c = results.counters;
    CHECK(c.mustBeZeroTotal() == 0);
    CHECK(c.np_model_ambiguous == 0);
    const auto& mix = results.model_mix;
    CHECK(mix.classified == kEvents);
    CHECK(results.np_hp_events + mix.observed_non_hp == kEvents);
    // About half of a flat window goes to CHIPS.
    CHECK_THAT(mix.expected_non_hp, WithinRel(0.5 * kEvents, 0.02));
    REQUIRE(mix.skipped.empty());
    CHECK(std::abs(mix.pull) < 4.0);
    REQUIRE(mix.chi2.has_value());
    CHECK(mix.ndf == static_cast<int>(checks::kModelMixBins));
    CHECK(mix.p_value > 1e-4);
    // Both oracles hold on their own sample.
    CHECK(mix.Ep_residuals.samples == mix.observed_non_hp);
    CHECK(mix.Ep_residuals.max_abs < 1e-9 * checks::kHpElasticMaxEnergy);
    CHECK(results.Ep_residuals.samples == results.np_hp_events);
    CHECK(results.Ep_residuals.max_abs < 1e-9 * checks::kHpElasticMaxEnergy);
    // HP events keep their offset of about m_e/m_p from the physical
    // kinematics (its exact range depends on En, about 5.33e-4 near 20 MeV);
    // a CHIPS event leaking into this sample would pull the minimum to zero.
    CHECK(results.Ep_physical_offset.relative_min > 5.0e-4);
    REQUIRE(results.legendre.moments.has_value());
    CHECK(results.legendre.moments->samples == results.np_hp_events);

    // A step at the centre of the window instead of the ramp has the same
    // mean share but the wrong shape: the per-bin chi^2 rejects it.
    const auto step = windowRun([](double En) { return En > 19.75 ? 1.0 : 0.0; }, kEvents);
    REQUIRE(step.model_mix.chi2.has_value());
    CHECK(step.model_mix.p_value < 1e-4);
    // No CHIPS events at all: the total is off by many standard deviations.
    const auto none = windowRun([](double) { return 0.0; }, kEvents);
    CHECK(none.model_mix.pull < -10.0);
}

TEST_CASE("a non-HP event below the window is a violation", "[check_accumulator]") {
    const double target_mass = kAwrH1 * kin::kNeutronMass;
    checks::CheckAccumulator accumulator(cubeSetup(target_mass), {model("lib", kCoefficients)});
    auto r = npRecord(kEn, 0.5 * kin::maxProtonEnergy(kEn), 30.0);
    r.Ep = kin::protonEnergy(kEn, r.theta_p);
    accumulator.add(r);
    const auto results = accumulator.finish();
    CHECK(results.counters.np_non_hp_outside_window == 1);
    CHECK(results.counters.mustBeZeroTotal() == 1);
    CHECK(results.np_hp_events == 0);
    // At 14 MeV no CHIPS events are expected: nothing to compare.
    CHECK_FALSE(results.model_mix.skipped.empty());
    CHECK_FALSE(results.model_mix.chi2.has_value());
    CHECK_FALSE(results.model_mix.chi2_skipped.empty());
    CHECK(results.model_mix.edges.empty());
}

TEST_CASE("checks with too few events are skipped, not failed", "[check_accumulator]") {
    checks::CheckAccumulator accumulator(cubeSetup(), {model("lib", kCoefficients)});
    const double Ep_max = kin::maxProtonEnergy(kEn);
    accumulator.add(npRecord(kEn, 0.2 * Ep_max, 10.0));
    accumulator.add(npRecord(kEn, 0.7 * Ep_max, 60.0));
    accumulator.add(npRecord(kEn, 0.9 * Ep_max, 80.0));
    const auto results = accumulator.finish();

    const auto& u = results.spectrum_u.per_library.at(0);
    CHECK_FALSE(u.comparison.has_value());
    CHECK_FALSE(u.skipped.empty());
    CHECK_FALSE(results.depth_z.per_library.at(0).comparison.has_value());
    const auto& f = results.fractions.per_library.at(0);
    CHECK_FALSE(f.comparison.has_value());
    CHECK_FALSE(f.skipped.empty());
    // Three spread cosines are enough for the moments themselves.
    CHECK(results.legendre.moments.has_value());
}

TEST_CASE("nuclide from HP cross sections", "[check_accumulator][hp]") {
    // Round numbers on overlapping grids, not physical values.
    const checks::HpCrossSection elastic{
        {1.0, 4.0, 4.0, 10.0},
        {0.6 * units::barn, 0.6 * units::barn, 0.7 * units::barn, 0.4 * units::barn}};
    const checks::HpCrossSection inelastic{{2.0, 6.0, 12.0},
                                           {0.0, 0.2 * units::barn, 0.2 * units::barn}};
    const checks::HpCrossSection capture{{0.5, 8.0}, {2.0e-5 * units::barn, 1.0e-5 * units::barn}};
    constexpr double kDensity = 1.0e16; // atoms / mm^3
    const checks::Nuclide h2 = checks::makeHpNuclide(checks::Target::Deuterium, "H-2", kDensity,
                                                     elastic, inelastic, capture);

    SECTION("the total is the sum of the parts over the common interval") {
        CHECK(h2.target == checks::Target::Deuterium);
        CHECK_THAT(h2.number_density, WithinRel(kDensity, kTolerance));
        CHECK_THAT(h2.total.x.front(), WithinRel(2.0, kTolerance));
        CHECK_THAT(h2.total.x.back(), WithinRel(8.0, kTolerance));
        // Points between the union nodes too, where only linearity makes it exact.
        for (const double En : {2.0, 2.5, 3.9, 4.5, 5.0, 6.0, 7.3, 8.0}) {
            const double sum = checks::evaluate(elastic, En) + checks::evaluate(inelastic, En) +
                               checks::evaluate(capture, En);
            CHECK_THAT(endf::evaluate(h2.total, En), WithinRel(sum, kTolerance));
            CHECK_THAT(endf::evaluate(h2.elastic, En),
                       WithinRel(checks::evaluate(elastic, En), kTolerance));
            CHECK_THAT(endf::evaluate(h2.capture, En),
                       WithinRel(checks::evaluate(capture, En), kTolerance));
        }
    }
    SECTION("a jump of a part is a jump of the total") {
        const auto first =
            std::find_if(h2.total.x.begin(), h2.total.x.end(), [](double x) { return !(x < 4.0); });
        REQUIRE(std::next(first) != h2.total.x.end());
        const auto i = static_cast<std::size_t>(std::distance(h2.total.x.begin(), first));
        CHECK_THAT(h2.total.x.at(i + 1), WithinRel(4.0, kTolerance));
        const double others = checks::evaluate(inelastic, 4.0) + checks::evaluate(capture, 4.0);
        CHECK_THAT(h2.total.y.at(i), WithinRel(0.6 * units::barn + others, kTolerance));
        CHECK_THAT(h2.total.y.at(i + 1), WithinRel(0.7 * units::barn + others, kTolerance));
    }
    SECTION("disjoint or malformed tables are refused") {
        const checks::HpCrossSection late{{20.0, 30.0}, {0.1 * units::barn, 0.1 * units::barn}};
        CHECK_THROWS_AS(checks::makeHpNuclide(checks::Target::Deuterium, "H-2", kDensity, elastic,
                                              late, capture),
                        checks::HpDataError);
        const checks::HpCrossSection decreasing{
            {1.0, 3.0, 2.0}, {0.1 * units::barn, 0.1 * units::barn, 0.1 * units::barn}};
        CHECK_THROWS_AS(checks::makeHpNuclide(checks::Target::Deuterium, "H-2", kDensity,
                                              decreasing, inelastic, capture),
                        endf::EndfError);
    }
}
