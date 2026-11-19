#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <cstdint>
#include <stdexcept>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/deposit_checks.hpp"

namespace checks = riptide::sim::checks;
namespace units = riptide::units;
using Catch::Matchers::WithinAbs;
using riptide::schema::EventClass;
using riptide::schema::FirstInteraction;
using riptide::sim::DepositRecord;
using riptide::sim::TruthRecord;

namespace {

constexpr double kEnMax = 20.0 * units::MeV;
constexpr double kHalfEdge = 50.0 * units::mm;
constexpr double kKB = 0.155 * units::mm / units::MeV;
constexpr std::size_t kBins = 20; // 1 MeV bins over [0, 20] MeV
constexpr double kTolerance = 1e-12;

checks::DepositSetup setup(double kB = kKB) {
    return {{kHalfEdge, kHalfEdge, kHalfEdge}, kEnMax, kB, kBins};
}

std::int32_t code(EventClass c) {
    return static_cast<std::int32_t>(c);
}

// A single n-p scatter at the centre whose proton stops inside the block
// after depositing exactly its energy, quenched to 80 % light.
struct NpEvent {
    TruthRecord truth;
    DepositRecord deposit;
};

NpEvent containedProton(double Ep) {
    NpEvent e;
    e.truth.En = 14.0 * units::MeV;
    e.truth.first_interaction = FirstInteraction::NpElastic;
    e.truth.n_np_scatters = 1;
    e.truth.Ep = Ep;
    e.truth.valid_np = true;
    e.deposit.Edep = Ep;
    e.deposit.Evis = 0.8 * Ep;
    e.deposit.proton.seen = true;
    e.deposit.proton.end = {1.0, 2.0, 3.0};
    e.deposit.proton.length = 1.0;
    e.deposit.proton.Edep = Ep;
    e.deposit.proton.Evis = 0.8 * Ep;
    return e;
}

checks::DepositResults one(const NpEvent& e, std::int32_t stored = code(EventClass::SingleNp),
                           double kB = kKB) {
    checks::DepositAccumulator acc(setup(kB));
    acc.add(e.truth, e.deposit, stored);
    return acc.finish();
}

} // namespace

TEST_CASE("makeDepositSetup reads the box, the source and kB", "[deposit_checks]") {
    riptide::sim::RunConfig config;
    config.geometry.scint_size = {100.0, 60.0, 40.0};
    config.source.En_max = kEnMax;
    CHECK_THROWS_AS(checks::makeDepositSetup(config), std::invalid_argument);
    config.scintillator = riptide::sim::ScintillatorConfig{kKB, "D18"};
    const checks::DepositSetup s = checks::makeDepositSetup(config, kBins);
    CHECK_THAT(s.half_size.x, WithinAbs(50.0, kTolerance));
    CHECK_THAT(s.half_size.y, WithinAbs(30.0, kTolerance));
    CHECK_THAT(s.half_size.z, WithinAbs(20.0, kTolerance));
    CHECK_THAT(s.En_max, WithinAbs(kEnMax, kTolerance));
    CHECK_THAT(s.birks_kB, WithinAbs(kKB, kTolerance));
    CHECK(s.bins == kBins);
}

TEST_CASE("DepositAccumulator rejects an invalid setup", "[deposit_checks]") {
    checks::DepositSetup s = setup();
    s.bins = 0;
    CHECK_THROWS_AS(checks::DepositAccumulator(s), std::invalid_argument);
    s = setup(-kKB);
    CHECK_THROWS_AS(checks::DepositAccumulator(s), std::invalid_argument);
    s = setup();
    s.En_max = 0.0;
    CHECK_THROWS_AS(checks::DepositAccumulator(s), std::invalid_argument);
    s = setup();
    s.half_size.y = 0.0;
    CHECK_THROWS_AS(checks::DepositAccumulator(s), std::invalid_argument);
}

TEST_CASE("a consistent contained proton passes every check", "[deposit_checks]") {
    const checks::DepositResults r = one(containedProton(5.5 * units::MeV));
    CHECK(r.events == 1);
    CHECK(r.np_events == 1);
    CHECK(r.counters.mustBeZeroTotal() == 0);
    CHECK(r.counters.p_Edep_above_Ep == 0);
    CHECK(r.classes.at(static_cast<std::size_t>(EventClass::SingleNp)) == 1);
    CHECK(r.contained_balance.samples == 1);
    CHECK_THAT(r.contained_balance.mean, WithinAbs(0.0, kTolerance));
    CHECK(r.Edep_counts.at(5) == 1);
    CHECK(r.Evis_counts.at(4) == 1); // 4.4 MeV
    CHECK(r.protons_per_Ep.at(5) == 1);
    CHECK(r.protons_exiting == 0);
}

TEST_CASE("an event with no interaction and no deposit is consistent", "[deposit_checks]") {
    checks::DepositAccumulator acc(setup());
    acc.add(TruthRecord{}, DepositRecord{}, code(EventClass::NoInteraction));
    const checks::DepositResults r = acc.finish();
    CHECK(r.counters.mustBeZeroTotal() == 0);
    CHECK(r.np_events == 0);
    CHECK(r.classes.at(0) == 1);
    // Zero deposits are not histogrammed: the spectra hold interacting events.
    for (const std::uint64_t n : r.Edep_counts) {
        CHECK(n == 0);
    }
}

TEST_CASE("each must-be-zero counter fires on its own violation", "[deposit_checks]") {
    const double Ep = 5.5 * units::MeV;

    SECTION("visible above deposited") {
        NpEvent e = containedProton(Ep);
        e.deposit.Evis = e.deposit.Edep + 1e-6;
        CHECK(one(e).counters.Evis_above_Edep == 1);
    }
    SECTION("proton visible above proton deposited") {
        NpEvent e = containedProton(Ep);
        e.deposit.proton.Evis = Ep + 1e-6;
        CHECK(one(e).counters.p_Evis_above_p_Edep == 1);
    }
    SECTION("proton deposit above the event total") {
        NpEvent e = containedProton(Ep);
        e.deposit.Edep = Ep - 1e-6;
        CHECK(one(e).counters.p_Edep_above_det_Edep == 1);
    }
    SECTION("negative deposit") {
        NpEvent e = containedProton(Ep);
        e.deposit.Evis = -1e-6;
        CHECK(one(e).counters.negative_deposit == 1);
    }
    SECTION("proton deposit above Ep beyond the slack") {
        NpEvent e = containedProton(Ep);
        e.deposit.proton.Edep = Ep + 2.0 * checks::kNuclearStoppingSlack;
        e.deposit.Edep = e.deposit.proton.Edep;
        const checks::DepositResults r = one(e);
        CHECK(r.counters.p_Edep_above_Ep_slack == 1);
        CHECK(r.counters.p_Edep_above_Ep == 0);
    }
    SECTION("proton without a track") {
        NpEvent e = containedProton(Ep);
        e.deposit.proton.length = 0.0;
        CHECK(one(e).counters.proton_without_track == 1);
    }
    SECTION("deposit without any interaction") {
        DepositRecord d;
        d.Edep = 1.0;
        d.Evis = 0.5;
        checks::DepositAccumulator acc(setup());
        acc.add(TruthRecord{}, d, code(EventClass::NoInteraction));
        CHECK(acc.finish().counters.deposit_without_interaction == 1);
    }
    SECTION("stored class differs from the counters") {
        const checks::DepositResults r = one(containedProton(Ep), code(EventClass::MultipleNp));
        CHECK(r.counters.class_mismatch == 1);
        // The class counts follow the counters, not the stored column.
        CHECK(r.classes.at(static_cast<std::size_t>(EventClass::SingleNp)) == 1);
    }
    SECTION("first interaction not among the counted ones") {
        NpEvent e = containedProton(Ep);
        e.truth.first_interaction = FirstInteraction::NCElastic;
        CHECK(one(e).counters.first_interaction_uncounted == 1);
    }
    SECTION("exiting proton ending inside the block") {
        NpEvent e = containedProton(Ep);
        e.deposit.proton.exits = true;
        CHECK(one(e).counters.exit_not_on_surface == 1);
    }
    SECTION("proton ending outside the block") {
        NpEvent e = containedProton(Ep);
        e.deposit.proton.end = {0.0, 0.0, kHalfEdge + 1.0};
        CHECK(one(e).counters.end_outside_box == 1);
    }
    SECTION("kB = 0 with quenched light") {
        const checks::DepositResults r = one(containedProton(Ep), code(EventClass::SingleNp), 0.0);
        CHECK(r.counters.kB_zero_visible_mismatch == 1);
        CHECK_THAT(r.kB_zero_max_abs_difference, WithinAbs(0.2 * Ep, kTolerance));
    }
}

TEST_CASE("kB = 0 with equal totals passes the visible-energy oracle", "[deposit_checks]") {
    NpEvent e = containedProton(5.5 * units::MeV);
    e.deposit.Evis = e.deposit.Edep;
    e.deposit.proton.Evis = e.deposit.proton.Edep;
    const checks::DepositResults r = one(e, code(EventClass::SingleNp), 0.0);
    CHECK(r.counters.mustBeZeroTotal() == 0);
    CHECK_THAT(r.kB_zero_max_abs_difference, WithinAbs(0.0, 0.0));
}

TEST_CASE("known physics is counted as expected-rare, not as a failure", "[deposit_checks]") {
    const double Ep = 5.5 * units::MeV;
    SECTION("nuclear-stopping excess within the slack") {
        NpEvent e = containedProton(Ep);
        e.deposit.proton.Edep = Ep + 2e-4;
        e.deposit.Edep = e.deposit.proton.Edep;
        const checks::DepositResults r = one(e);
        CHECK(r.counters.mustBeZeroTotal() == 0);
        CHECK(r.counters.p_Edep_above_Ep == 1);
        CHECK(r.contained_balance.positive == 1);
        CHECK_THAT(r.contained_balance.max, WithinAbs(2e-4, kTolerance));
    }
    SECTION("contained proton stopping on the surface") {
        NpEvent e = containedProton(Ep);
        e.deposit.proton.end = {kHalfEdge, 0.0, 0.0};
        const checks::DepositResults r = one(e);
        CHECK(r.counters.mustBeZeroTotal() == 0);
        CHECK(r.counters.contained_end_on_surface == 1);
    }
}

TEST_CASE("exiting and interacting protons leave the energy balance", "[deposit_checks]") {
    const double Ep = 15.5 * units::MeV;
    checks::DepositAccumulator acc(setup());

    NpEvent exiting = containedProton(Ep);
    exiting.deposit.proton.exits = true;
    exiting.deposit.proton.end = {0.0, -kHalfEdge, 10.0};
    exiting.deposit.proton.Edep = 3.0;
    exiting.deposit.Edep = 3.0;
    exiting.deposit.proton.Evis = 2.0;
    exiting.deposit.Evis = 2.0;
    acc.add(exiting.truth, exiting.deposit, code(EventClass::SingleNp));

    NpEvent interacting = containedProton(Ep);
    interacting.deposit.proton.n_hadronic = 1;
    interacting.deposit.proton.Edep = 9.0;
    interacting.deposit.proton.Evis = 7.0;
    acc.add(interacting.truth, interacting.deposit, code(EventClass::SingleNp));

    const checks::DepositResults r = acc.finish();
    CHECK(r.counters.mustBeZeroTotal() == 0);
    CHECK(r.contained_balance.samples == 0);
    CHECK(r.protons_exiting == 1);
    CHECK(r.protons_hadronic == 1);
    CHECK(r.protons_per_Ep.at(15) == 2);
    CHECK(r.exiting_per_Ep.at(15) == 1);
}

TEST_CASE("deposits above En_max go to the overflow", "[deposit_checks]") {
    // A capture adds its Q-value, so the total can exceed the source energy.
    TruthRecord truth;
    truth.first_interaction = FirstInteraction::Other;
    truth.n_other = 1;
    DepositRecord d;
    d.Edep = kEnMax + 2.0;
    d.Evis = kEnMax; // the upper edge itself stays in the last bin
    checks::DepositAccumulator acc(setup());
    acc.add(truth, d, code(EventClass::Other));
    const checks::DepositResults r = acc.finish();
    CHECK(r.counters.mustBeZeroTotal() == 0);
    CHECK(r.Edep_overflow == 1);
    CHECK(r.Evis_overflow == 0);
    CHECK(r.Evis_counts.back() == 1);
}
