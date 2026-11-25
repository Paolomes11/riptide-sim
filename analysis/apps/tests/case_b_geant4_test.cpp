#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <stdexcept>

#include "riptide/apps/case_b_geant4.hpp"
#include "riptide/io/run_info.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;
using riptide::apps::CliError;
using riptide::apps::parseDirection;
using riptide::geometry::Direction;

constexpr double kUnitTolerance = 1e-15;

riptide::EventSample twoEvents() {
    const double s = 1.0 / std::sqrt(2.0);
    riptide::EventSample events;
    events.event_id = {7, 9};
    events.truth_En = {14.0, 15.0};
    events.truth_omega_n_x = {1.0, 1.0};
    events.truth_omega_n_y = {0.0, 0.0};
    events.truth_omega_n_z = {0.0, 0.0};
    events.truth_Ep = {5.01, 6.01};
    events.truth_Ep_phys = {5.0, 6.0};
    events.truth_theta_p = {0.5, 0.6};
    events.truth_p_dir_x = {s, 0.0};
    events.truth_p_dir_y = {s, 1.0};
    events.truth_p_dir_z = {0.0, 0.0};
    events.reco_Ep = {5.1, 5.9};
    events.reco_theta_p = {0.45, 0.65};
    events.reco_p_dir_x = {0.0, 0.0};
    events.reco_p_dir_y = {0.0, s};
    events.reco_p_dir_z = {1.0, s};
    return events;
}

nlohmann::json runRecord(const nlohmann::json& direction) {
    return {{"config", {{"source", {{"beam", "far_field"}, {"direction", direction}}}}}};
}

} // namespace

TEST_CASE("the Case B sample takes the physical energy and both tracks", "[case_b_geant4]") {
    const auto sample = riptide::apps::caseBSampleFromEvents(twoEvents());
    REQUIRE(sample.size() == 2);
    REQUIRE(sample.truth_track.size() == 2);
    REQUIRE(sample.reco_track.size() == 2);
    CHECK_THAT(sample.truth_En[1], WithinULP(15.0, 0));
    CHECK_THAT(sample.truth_Ep[0], WithinULP(5.0, 0));
    CHECK_THAT(sample.truth_theta_p[1], WithinULP(0.6, 0));
    CHECK_THAT(sample.reco_Ep[0], WithinULP(5.1, 0));
    CHECK_THAT(sample.truth_track[1].y, WithinULP(1.0, 0));
    CHECK_THAT(sample.reco_track[0].z, WithinULP(1.0, 0));
    CHECK_THAT(sample.reco_track[1].y, WithinULP(1.0 / std::sqrt(2.0), 0));
}

TEST_CASE("a partly filled or trackless Case B event sample is refused", "[case_b_geant4]") {
    auto events = twoEvents();
    events.reco_p_dir_z.pop_back();
    CHECK_THROWS_AS(riptide::apps::caseBSampleFromEvents(events), std::invalid_argument);
    events = twoEvents();
    events.truth_Ep_phys.clear();
    CHECK_THROWS_AS(riptide::apps::caseBSampleFromEvents(events), std::invalid_argument);
    // The zero vector marks an event without a proton track.
    events = twoEvents();
    events.truth_p_dir_y[1] = 0.0;
    CHECK_THROWS_AS(riptide::apps::caseBSampleFromEvents(events), std::invalid_argument);
    events = twoEvents();
    events.reco_p_dir_z[0] = 1.0 + 1e-6;
    CHECK_THROWS_AS(riptide::apps::caseBSampleFromEvents(events), std::invalid_argument);
}

TEST_CASE("the direction spec gives a normalised unit vector", "[case_b_geant4]") {
    const Direction x = parseDirection("1,0,0");
    CHECK_THAT(x.x, WithinULP(1.0, 0));
    CHECK_THAT(x.y, WithinULP(0.0, 0));
    // The oblique direction of D45, (3, 2, 1) / sqrt(14), as written in the config.
    const Direction obl =
        parseDirection("0.8017837257372732,0.5345224838248488,0.2672612419124244");
    CHECK_THAT(obl.x, WithinAbs(3.0 / std::sqrt(14.0), kUnitTolerance));
    CHECK_THAT(obl.y, WithinAbs(2.0 / std::sqrt(14.0), kUnitTolerance));
    CHECK_THAT(obl.z, WithinAbs(1.0 / std::sqrt(14.0), kUnitTolerance));
    // Within the norm tolerance the vector is normalised.
    const Direction z = parseDirection("0,0,1.0000005");
    CHECK_THAT(z.z, WithinULP(1.0, 0));
}

TEST_CASE("the direction spec refuses malformed specs", "[case_b_geant4]") {
    for (const char* spec : {"", "1", "1,0", "1,0,0,0", "1,,0", ",0,1", "1,0,", "1;0;0", "x,0,0",
                             "0,0,0", "2,0,0", "0,0,1.00001", "nan,0,0", "inf,0,0", "1,0,0 "}) {
        INFO(spec);
        CHECK_THROWS_AS(parseDirection(spec), CliError);
    }
}

TEST_CASE("the fixed direction of a run comes from its config", "[case_b_geant4]") {
    const Direction z = riptide::apps::fixedDirectionFromRunRecord(runRecord({0.0, 0.0, 1.0}));
    CHECK_THAT(z.z, WithinULP(1.0, 0));
    const Direction obl = riptide::apps::fixedDirectionFromRunRecord(
        runRecord({0.8017837257372732, 0.5345224838248488, 0.2672612419124244}));
    CHECK_THAT(obl.x, WithinAbs(3.0 / std::sqrt(14.0), kUnitTolerance));
}

TEST_CASE("a run without a fixed direction is refused", "[case_b_geant4]") {
    using riptide::apps::fixedDirectionFromRunRecord;
    using riptide::io::IoError;
    CHECK_THROWS_AS(fixedDirectionFromRunRecord(runRecord("isotropic")), IoError);
    CHECK_THROWS_AS(fixedDirectionFromRunRecord(runRecord({1.0, 0.0})), IoError);
    CHECK_THROWS_AS(fixedDirectionFromRunRecord(runRecord({1.0, 0.0, "0"})), IoError);
    CHECK_THROWS_AS(fixedDirectionFromRunRecord(runRecord({2.0, 0.0, 0.0})), IoError);
    CHECK_THROWS_AS(fixedDirectionFromRunRecord(nlohmann::json{{"config", {{"source", {}}}}}),
                    IoError);
    CHECK_THROWS_AS(fixedDirectionFromRunRecord(nlohmann::json::object()), IoError);
}

TEST_CASE("every event must carry the configured direction", "[case_b_geant4]") {
    auto events = twoEvents();
    CHECK_NOTHROW(riptide::apps::requireFixedDirection(events, {1.0, 0.0, 0.0}));
    CHECK_THROWS_AS(riptide::apps::requireFixedDirection(events, {0.0, 0.0, 1.0}),
                    std::invalid_argument);
    events.truth_omega_n_y[1] = 1e-6;
    CHECK_THROWS_AS(riptide::apps::requireFixedDirection(events, {1.0, 0.0, 0.0}),
                    std::invalid_argument);
    events = twoEvents();
    events.truth_omega_n_z.pop_back();
    CHECK_THROWS_AS(riptide::apps::requireFixedDirection(events, {1.0, 0.0, 0.0}),
                    std::invalid_argument);
}
