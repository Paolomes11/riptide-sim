#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/check_accumulator.hpp"
#include "riptide/sim/expected_distributions.hpp"
#include "riptide/sim/run_prediction.hpp"

namespace checks = riptide::sim::checks;
namespace endf = riptide::sim::endf;
namespace kin = riptide::kinematics;
namespace units = riptide::units;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using riptide::schema::FirstInteraction;

namespace {

constexpr double kEn = 14.0 * units::MeV;
constexpr double kTolerance = 1e-12;
constexpr std::uint64_t kEvents = 10'000'000;

// Plausible plastic-scintillator medium, energy independent so that every
// spectrum average equals its single-energy value.
constexpr double kHydrogenDensity = 5.0e19; // atoms / mm^3
constexpr double kCarbonDensity = 4.0e19;   // atoms / mm^3
constexpr double kHydrogenElastic = 0.7 * units::barn;
constexpr double kCarbonElastic = 0.8 * units::barn;
constexpr double kCarbonCapture = 1.0e-5 * units::barn;
constexpr double kCarbonTotal = 1.3 * units::barn;

// Target mass of H-1 in the G4NDL 4.7.1 elastic file, in neutron masses.
constexpr double kAwrH1 = 0.999167;

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

// Pencil beam along +z through a 100 mm cube: L = 100 mm.
checks::CheckSetup cubeSetup(double En_min, double En_max,
                             double target_mass = kin::kProtonMass) {
    riptide::sim::RunConfig config;
    config.geometry.scint_size = {100.0, 100.0, 100.0};
    config.source.En_min = En_min;
    config.source.En_max = En_max;
    config.source.direction = {0.0, 0.0, 1.0};
    config.source.position = {0.0, 0.0, -100.0};
    return checks::makeCheckSetup(config, target_mass, checks::kDefaultBins,
                                  checks::kDefaultBins);
}

std::size_t index(FirstInteraction c) {
    return static_cast<std::size_t>(c);
}

} // namespace

TEST_CASE("gaussLegendre integrates polynomials up to degree 2n - 1", "[run_prediction]") {
    for (std::size_t n = 1; n <= 12; ++n) {
        const auto rule = checks::gaussLegendre(n);
        REQUIRE(rule.nodes.size() == n);
        for (std::size_t i = 1; i < n; ++i) {
            CHECK(rule.nodes.at(i) > rule.nodes.at(i - 1));
        }
        for (std::size_t k = 0; k <= 2 * n - 1; ++k) {
            double integral = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                integral += rule.weights.at(i) * std::pow(rule.nodes.at(i), static_cast<double>(k));
            }
            const double exact = k % 2 == 0 ? 2.0 / static_cast<double>(k + 1) : 0.0;
            CHECK_THAT(integral, WithinAbs(exact, kTolerance));
        }
    }
    CHECK_THROWS_AS(checks::gaussLegendre(0), std::invalid_argument);
}

TEST_CASE("sourceQuadrature: one node for a line, midpoints for a flat source",
          "[run_prediction]") {
    const auto mono = checks::sourceQuadrature(kEn, kEn, 8);
    REQUIRE(mono.nodes.size() == 1);
    CHECK_THAT(mono.nodes.front(), WithinRel(kEn, kTolerance));
    CHECK_THAT(mono.weights.front(), WithinRel(1.0, kTolerance));

    const auto flat = checks::sourceQuadrature(10.0, 20.0, 4);
    REQUIRE(flat.nodes.size() == 4);
    const std::vector<double> midpoints{11.25, 13.75, 16.25, 18.75};
    double total = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK_THAT(flat.nodes.at(i), WithinRel(midpoints.at(i), kTolerance));
        total += flat.weights.at(i);
    }
    CHECK_THAT(total, WithinRel(1.0, kTolerance));

    CHECK_THROWS_AS(checks::sourceQuadrature(20.0, 10.0, 4), std::invalid_argument);
    CHECK_THROWS_AS(checks::sourceQuadrature(0.0, 10.0, 4), std::invalid_argument);
    CHECK_THROWS_AS(checks::sourceQuadrature(10.0, 20.0, 0), std::invalid_argument);
}

TEST_CASE("legendreSecondMoment matches closed forms", "[run_prediction]") {
    // Isotropic: <P_l^2> = 1 / (2l + 1).
    for (int l = 1; l <= 8; ++l) {
        CHECK_THAT(checks::legendreSecondMoment({}, l),
                   WithinRel(1.0 / (2.0 * l + 1.0), kTolerance));
    }
    // mu^2 = (2 P_2 + 1) / 3, so <P_1^2> = (2 a_2 + 1) / 3 for any f.
    const std::vector<double> a{0.1, 0.05, -0.02};
    CHECK_THAT(checks::legendreSecondMoment(a, 1), WithinRel((2.0 * a.at(1) + 1.0) / 3.0, 1e-12));
    CHECK_THROWS_AS(checks::legendreSecondMoment(a, 0), std::invalid_argument);
}

TEST_CASE("normalTwoSidedQuantile inverts the two-sided normal tail", "[run_prediction]") {
    // P(|Z| > 3) and P(|Z| > 1.959963984540054) = 0.05.
    CHECK_THAT(checks::normalTwoSidedQuantile(std::erfc(3.0 / std::sqrt(2.0))),
               WithinRel(3.0, 1e-12));
    CHECK_THAT(checks::normalTwoSidedQuantile(0.05), WithinRel(1.959963984540054, 1e-12));
    CHECK_THROWS_AS(checks::normalTwoSidedQuantile(0.0), std::invalid_argument);
    CHECK_THROWS_AS(checks::normalTwoSidedQuantile(1.0), std::invalid_argument);
}

TEST_CASE("predictRun of a monoenergetic source equals the single-energy formulas",
          "[run_prediction]") {
    const auto setup = cubeSetup(kEn, kEn);
    const auto m = model("lib", kCoefficients);
    const auto run = checks::predictRun(setup, {m}, kEvents, 8);
    CHECK(run.energy_nodes == 1);
    CHECK_FALSE(run.overlaps_model_window);
    CHECK_THAT(run.L, WithinRel(100.0, kTolerance));
    REQUIRE(run.libraries.size() == 1);
    const auto& lib = run.libraries.front();

    const auto first = checks::firstInteraction(m.mediumAt(kEn), setup.chord.length());
    CHECK_THAT(lib.sigma_total, WithinRel(first.sigma_total, kTolerance));
    for (std::size_t c = 0; c < checks::kFirstInteractionCategories; ++c) {
        CHECK_THAT(lib.probabilities.at(c), WithinAbs(first.unconditional.at(c), kTolerance));
    }
    CHECK_THAT(lib.hp_np_probability,
               WithinRel(first.unconditional.at(index(FirstInteraction::NpElastic)), kTolerance));
    CHECK_THAT(lib.non_hp_np_probability, WithinAbs(0.0, kTolerance));

    REQUIRE(lib.legendre.size() == kCoefficients.size());
    for (std::size_t l = 0; l < kCoefficients.size(); ++l) {
        CHECK_THAT(lib.legendre.at(l), WithinRel(kCoefficients.at(l), kTolerance));
        const double second = checks::legendreSecondMoment(kCoefficients, static_cast<int>(l + 1));
        CHECK_THAT(lib.legendre_spread.at(l),
                   WithinRel(std::sqrt(second - kCoefficients.at(l) * kCoefficients.at(l)),
                             kTolerance));
    }
    const auto u = checks::protonEnergyBinProbabilities(kCoefficients,
                                                        [] {
                                                            std::vector<double> edges;
                                                            for (std::size_t i = 0;
                                                                 i <= checks::kDefaultBins; ++i) {
                                                                edges.push_back(
                                                                    static_cast<double>(i) /
                                                                    checks::kDefaultBins);
                                                            }
                                                            return edges;
                                                        }());
    REQUIRE(lib.u_probabilities.size() == u.size());
    for (std::size_t b = 0; b < u.size(); ++b) {
        CHECK_THAT(lib.u_probabilities.at(b), WithinAbs(u.at(b), kTolerance));
    }

    // With the bare proton as HP target the HP kinematics is the physical one.
    CHECK_THAT(run.offset_relative_min, WithinAbs(0.0, kTolerance));
    CHECK_THAT(run.offset_relative_max, WithinAbs(0.0, kTolerance));

    const auto errors = checks::expectedErrors(run, 0);
    const double hp_events = static_cast<double>(kEvents) * lib.hp_np_probability;
    CHECK_THAT(errors.hp_np_events, WithinRel(hp_events, kTolerance));
    CHECK_THAT(errors.legendre_sigma.at(0),
               WithinRel(lib.legendre_spread.at(0) / std::sqrt(hp_events), kTolerance));
    const double p_np = lib.probabilities.at(index(FirstInteraction::NpElastic));
    CHECK_THAT(errors.fraction_sigma.at(index(FirstInteraction::NpElastic)),
               WithinRel(std::sqrt(p_np * (1.0 - p_np) / static_cast<double>(kEvents)),
                         kTolerance));
}

TEST_CASE("predictRun: the H-1 atomic mass shifts Ep upward by a few 1e-4",
          "[run_prediction]") {
    const auto setup = cubeSetup(kEn, kEn, kAwrH1 * kin::kNeutronMass);
    const auto run = checks::predictRun(setup, {model("lib", kCoefficients)}, kEvents, 1);
    // A lighter target takes more energy at a given lab angle; the order is
    // the relative mass difference, a few parts in 1e4.
    CHECK(run.offset_relative_min > 0.0);
    CHECK(run.offset_relative_max < 1.0e-3);
    CHECK(run.offset_relative_max >= run.offset_relative_min);
}

TEST_CASE("predictRun of a flat source: constant data average to themselves",
          "[run_prediction]") {
    const auto mono = checks::predictRun(cubeSetup(kEn, kEn), {model("lib", kCoefficients)},
                                         kEvents, 1);
    const auto flat = checks::predictRun(cubeSetup(10.0, 19.0), {model("lib", kCoefficients)},
                                         kEvents, 16);
    CHECK(flat.energy_nodes == 16);
    CHECK_FALSE(flat.overlaps_model_window);
    const auto& a = mono.libraries.front();
    const auto& b = flat.libraries.front();
    for (std::size_t c = 0; c < checks::kFirstInteractionCategories; ++c) {
        CHECK_THAT(b.probabilities.at(c), WithinAbs(a.probabilities.at(c), kTolerance));
    }
    const auto shift = checks::predictionShift(a, b, checks::expectedErrors(mono, 0), kEvents);
    CHECK_THAT(shift.fractions, WithinAbs(0.0, 1e-6));
    CHECK_THAT(shift.legendre, WithinAbs(0.0, 1e-6));
    CHECK_THAT(shift.u_noncentrality, WithinAbs(0.0, 1e-9));
}

TEST_CASE("predictRun across the overlap window splits n-p between HP and CHIPS",
          "[run_prediction]") {
    const auto run = checks::predictRun(cubeSetup(checks::kChipsElasticMinEnergy,
                                                  checks::kHpElasticMaxEnergy),
                                        {model("lib", kCoefficients)}, kEvents, 10);
    CHECK(run.overlaps_model_window);
    const auto& lib = run.libraries.front();
    const double np = lib.probabilities.at(index(FirstInteraction::NpElastic));
    // The CHIPS share is a linear ramp over the window: its mean is 1/2, and
    // the midpoint rule is exact for it.
    CHECK_THAT(lib.non_hp_np_probability, WithinRel(0.5 * np, 1e-12));
    CHECK_THAT(lib.hp_np_probability, WithinRel(0.5 * np, 1e-12));

    // At the upper edge every n-p event is CHIPS: nothing to predict for HP.
    CHECK_THROWS_AS(checks::predictRun(cubeSetup(checks::kHpElasticMaxEnergy,
                                                 checks::kHpElasticMaxEnergy),
                                       {model("lib", kCoefficients)}, kEvents, 1),
                    std::invalid_argument);
}

TEST_CASE("predictRun pads shorter Legendre expansions with zeros", "[run_prediction]") {
    const auto run = checks::predictRun(cubeSetup(kEn, kEn),
                                        {model("short", kCoefficients),
                                         model("long", {0.1, 0.05, 0.01})},
                                        kEvents, 1);
    REQUIRE(run.libraries.size() == 2);
    REQUIRE(run.libraries.at(0).legendre.size() == 3);
    CHECK_THAT(run.libraries.at(0).legendre.at(2), WithinAbs(0.0, kTolerance));
    CHECK_THAT(run.libraries.at(1).legendre.at(2), WithinRel(0.01, kTolerance));

    const auto errors = checks::expectedErrors(run, 0);
    const auto shift =
        checks::predictionShift(run.libraries.at(0), run.libraries.at(1), errors, kEvents);
    CHECK_THAT(shift.legendre, WithinRel(0.01 / errors.legendre_sigma.at(2), 1e-12));
    CHECK_THAT(shift.fractions, WithinAbs(0.0, kTolerance));
    CHECK(shift.u_noncentrality > 0.0);

    CHECK_THROWS_AS(checks::predictRun(cubeSetup(kEn, kEn), {}, kEvents, 1),
                    std::invalid_argument);
    CHECK_THROWS_AS(checks::predictRun(cubeSetup(kEn, kEn), {model("lib", kCoefficients)}, 0, 1),
                    std::invalid_argument);
}

TEST_CASE("testPlan counts the tests and splits alpha", "[run_prediction]") {
    constexpr double kAlpha = 0.01;
    const auto mono =
        checks::predictRun(cubeSetup(kEn, kEn), {model("lib", kCoefficients)}, kEvents, 1);
    const auto plan = checks::testPlan(mono, kAlpha);
    CHECK(plan.chi2_tests == 3);
    CHECK(plan.pull_tests == 2);
    CHECK_THAT(plan.alpha_per_test, WithinRel(kAlpha / 5.0, kTolerance));
    CHECK_THAT(std::erfc(plan.pull_threshold / std::sqrt(2.0)),
               WithinRel(plan.alpha_per_test, 1e-9));

    const auto edge = checks::predictRun(
        cubeSetup(checks::kChipsElasticMinEnergy, checks::kHpElasticMaxEnergy),
        {model("lib", kCoefficients)}, kEvents, 10);
    const auto edge_plan = checks::testPlan(edge, kAlpha);
    CHECK(edge_plan.chi2_tests == 4);
    CHECK(edge_plan.pull_tests == 3);

    CHECK_THROWS_AS(checks::testPlan(mono, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(checks::testPlan(mono, 1.0), std::invalid_argument);
}

TEST_CASE("predictionJson carries the run, the plan and the library shifts",
          "[run_prediction]") {
    const auto run = checks::predictRun(cubeSetup(kEn, kEn),
                                        {model("a", kCoefficients), model("b", {0.11, 0.05})},
                                        kEvents, 1);
    const auto plan = checks::testPlan(run, 0.01);
    const auto doc = checks::predictionJson(run, checks::PredictionShift{}, plan);
    CHECK(doc.at("events").get<std::uint64_t>() == kEvents);
    CHECK(doc.at("test_plan").at("pull_tests").get<int>() == 2);
    REQUIRE(doc.at("libraries").size() == 2);
    CHECK_FALSE(doc.at("libraries").at(0).contains("shift_from_first_in_sigma"));
    CHECK(doc.at("libraries").at(1).at("shift_from_first_in_sigma").at("legendre").get<double>() >
          0.0);
    CHECK(doc.at("libraries").at(0).at("probabilities").contains("np_elastic"));
}
