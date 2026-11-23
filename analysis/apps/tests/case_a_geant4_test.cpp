#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <stdexcept>

#include "riptide/apps/case_a_geant4.hpp"

namespace {

using riptide::apps::CliError;
using riptide::apps::parseTruth;
using riptide::apps::rowExperiments;
using riptide::io::SpectrumMode;

riptide::EventSample twoEvents() {
    riptide::EventSample events;
    events.event_id = {7, 9};
    events.truth_En = {14.0, 15.0};
    events.truth_Ep = {5.01, 6.01};
    events.truth_Ep_phys = {5.0, 6.0};
    events.truth_theta_p = {0.5, 0.6};
    events.reco_Ep = {5.1, 5.9};
    events.reco_theta_p = {0.45, 0.65};
    return events;
}

} // namespace

TEST_CASE("the Case A sample takes the physical proton energy", "[case_a_geant4]") {
    const auto sample = riptide::apps::pseudoSampleFromEvents(twoEvents());
    REQUIRE(sample.size() == 2);
    CHECK_THAT(sample.truth_En[1], Catch::Matchers::WithinULP(15.0, 0));
    CHECK_THAT(sample.truth_Ep[0], Catch::Matchers::WithinULP(5.0, 0));
    CHECK_THAT(sample.truth_theta_p[1], Catch::Matchers::WithinULP(0.6, 0));
    CHECK_THAT(sample.reco_Ep[0], Catch::Matchers::WithinULP(5.1, 0));
    CHECK_THAT(sample.reco_theta_p[1], Catch::Matchers::WithinULP(0.65, 0));
}

TEST_CASE("a partly filled event sample is refused", "[case_a_geant4]") {
    auto events = twoEvents();
    events.reco_Ep.pop_back();
    CHECK_THROWS_AS(riptide::apps::pseudoSampleFromEvents(events), std::invalid_argument);
    events = twoEvents();
    events.truth_Ep_phys.clear();
    CHECK_THROWS_AS(riptide::apps::pseudoSampleFromEvents(events), std::invalid_argument);
}

TEST_CASE("the truth spec gives a mono or a flat spectrum", "[case_a_geant4]") {
    const auto mono = parseTruth("mono:14");
    CHECK(mono.mode == SpectrumMode::Mono);
    CHECK_THAT(mono.En_MeV, Catch::Matchers::WithinULP(14.0, 0));
    const auto flat = parseTruth("flat:10:19.5");
    CHECK(flat.mode == SpectrumMode::Flat);
    CHECK_THAT(flat.En_min_MeV, Catch::Matchers::WithinULP(10.0, 0));
    CHECK_THAT(flat.En_max_MeV, Catch::Matchers::WithinULP(19.5, 0));
}

TEST_CASE("the truth spec refuses malformed specs", "[case_a_geant4]") {
    for (const char* spec :
         {"", "mono", "mono:", "mono:-1", "mono:0", "mono:inf", "mono:nan", "mono:14:15",
          "mono:14x", "flat:10", "flat:20:10", "flat:10:10", "flat:10:20:30", "x:1", "MONO:14"}) {
        INFO(spec);
        CHECK_THROWS_AS(parseTruth(spec), CliError);
    }
}

TEST_CASE("a row has the planned experiments or as many as fit", "[case_a_geant4]") {
    // Planned M when the pool is large enough.
    CHECK(rowExperiments(400, 1000000, 10, false) == 400);
    // Fewer when not: 1000 events give 100 experiments of 10.
    CHECK(rowExperiments(400, 1000, 10, false) == 100);
    // The pilot scales M, at least two.
    CHECK(rowExperiments(400, 1000000, 10, true) == 40);
    CHECK(rowExperiments(15, 1000000, 10, true) == 2);
    // Fewer than two disjoint experiments: the row is skipped.
    CHECK(rowExperiments(40, 1999, 1000, false) == 0);
    CHECK(rowExperiments(40, 2000, 1000, false) == 2);
    CHECK(rowExperiments(40, 999, 1000, true) == 0);
    CHECK_THROWS_AS(rowExperiments(40, 1000, 0, false), std::invalid_argument);
}

TEST_CASE("the pilot divides counts by ten, at least two", "[case_a_geant4]") {
    CHECK(riptide::apps::pilotScaled(20000, true) == 2000);
    CHECK(riptide::apps::pilotScaled(20000, false) == 20000);
    CHECK(riptide::apps::pilotScaled(5, true) == 2);
}

TEST_CASE("a mono sample must carry one energy", "[case_a_geant4]") {
    riptide::stats::PseudoSample sample;
    sample.truth_En = {14.0, 14.0, 14.0};
    CHECK_NOTHROW(riptide::apps::requireMonoTruth(sample, 14.0));
    CHECK_THROWS_AS(riptide::apps::requireMonoTruth(sample, 15.0), std::invalid_argument);
    sample.truth_En.push_back(14.001);
    CHECK_THROWS_AS(riptide::apps::requireMonoTruth(sample, 14.0), std::invalid_argument);
}
