#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstdint>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/io/acceptance_record.hpp"
#include "riptide/io/selection.hpp"

namespace {

namespace io = riptide::io;
using Catch::Matchers::WithinULP;
using riptide::io::IoError;

constexpr double kHalfPi = 0.5 * riptide::units::pi;

riptide::EventSample sampleOf(const std::vector<double>& En, const std::vector<double>& theta) {
    riptide::EventSample s;
    for (std::size_t i = 0; i < En.size(); ++i) {
        s.event_id.push_back(static_cast<std::int32_t>(i));
        s.truth_En.push_back(En[i]);
        s.truth_theta_p.push_back(theta[i]);
    }
    return s;
}

} // namespace

TEST_CASE("stepped edges end at hi, with a shorter last bin", "[acceptance_record]") {
    const std::vector<double> edges = io::steppedEdges(10.0, 19.5, io::kAcceptanceEnStep);
    REQUIRE(edges.size() == 11);
    CHECK_THAT(edges.front(), WithinULP(10.0, 0));
    CHECK_THAT(edges[9], WithinULP(19.0, 0));
    CHECK_THAT(edges.back(), WithinULP(19.5, 0));
    // An exact multiple has no sliver at the end.
    const std::vector<double> exact = io::steppedEdges(10.0, 20.0, 1.0);
    REQUIRE(exact.size() == 11);
    CHECK_THAT(exact.back(), WithinULP(20.0, 0));
    CHECK_THROWS_AS(io::steppedEdges(1.0, 1.0, 1.0), IoError);
    CHECK_THROWS_AS(io::steppedEdges(1.0, 2.0, 0.0), IoError);
}

TEST_CASE("uniform edges cover [lo, hi] exactly", "[acceptance_record]") {
    const std::vector<double> edges = io::uniformEdges(0.0, kHalfPi, io::kAcceptanceThetaBins);
    REQUIRE(edges.size() == io::kAcceptanceThetaBins + 1);
    CHECK_THAT(edges.back(), WithinULP(kHalfPi, 0));
    CHECK_THAT(edges[1], WithinULP(kHalfPi / 18.0, 1));
    CHECK_THROWS_AS(io::uniformEdges(0.0, 1.0, 0), IoError);
}

TEST_CASE("the table counts each event once, the last bin closed", "[acceptance_record]") {
    const std::vector<double> En_edges{10.0, 11.0, 12.0};
    const std::vector<double> theta_edges{0.0, 0.5, kHalfPi};
    // Edges: lower edge in, interior edge to the upper bin, top edges closed.
    const auto sample = sampleOf({10.0, 11.0, 12.0, 10.5, 11.5}, {0.0, 0.5, kHalfPi, 0.2, 1.0});
    const io::AcceptanceTable t =
        io::fillAcceptanceTable(sample, En_edges, theta_edges, io::EventHalf::All);
    REQUIRE(t.counts.size() == 4);
    CHECK_THAT(t.counts[0], WithinULP(2.0, 0)); // (10, 0) and (10.5, 0.2)
    CHECK_THAT(t.counts[1], WithinULP(0.0, 0));
    CHECK_THAT(t.counts[2], WithinULP(0.0, 0));
    CHECK_THAT(t.counts[3], WithinULP(3.0, 0)); // (11, 0.5), (12, pi/2), (11.5, 1)

    const io::AcceptanceTable even =
        io::fillAcceptanceTable(sample, En_edges, theta_edges, io::EventHalf::Even);
    const io::AcceptanceTable odd =
        io::fillAcceptanceTable(sample, En_edges, theta_edges, io::EventHalf::Odd);
    for (std::size_t i = 0; i < t.counts.size(); ++i) {
        CHECK_THAT(even.counts[i] + odd.counts[i], WithinULP(t.counts[i], 0));
    }
    CHECK_THAT(even.counts[3], WithinULP(2.0, 0)); // events 2 and 4

    CHECK_THROWS_AS(
        io::fillAcceptanceTable(sampleOf({9.99}, {0.1}), En_edges, theta_edges, io::EventHalf::All),
        IoError);
    CHECK_THROWS_AS(
        io::fillAcceptanceTable(sampleOf({10.5}, {1.6}), En_edges, theta_edges, io::EventHalf::All),
        IoError);
    CHECK_THROWS_AS(io::fillAcceptanceTable(sample, {10.0}, theta_edges, io::EventHalf::All),
                    IoError);
    CHECK_THROWS_AS(
        io::fillAcceptanceTable(sample, {10.0, 10.0, 12.0}, theta_edges, io::EventHalf::All),
        IoError);
}

TEST_CASE("event halves parse by name", "[acceptance_record]") {
    CHECK(io::parseEventHalf("all") == io::EventHalf::All);
    CHECK(io::parseEventHalf("even") == io::EventHalf::Even);
    CHECK(io::parseEventHalf("odd") == io::EventHalf::Odd);
    CHECK_THROWS_AS(io::parseEventHalf("half"), IoError);
}

TEST_CASE("acceptance.json round-trips its table", "[acceptance_record]") {
    const auto sample = sampleOf({10.2, 11.7, 11.9}, {0.1, 1.2, 0.3});
    const io::AcceptanceTable t = io::fillAcceptanceTable(sample, {10.0, 11.0, 12.0},
                                                          {0.0, 0.5, kHalfPi}, io::EventHalf::All);
    io::RunInfo run;
    run.run_id = "run_a";
    run.seed = 20260908;
    const io::Selection selection = io::Selection::parse("single_np,has_proton,contained");
    const io::CutFlow flow(selection.cuts().size());
    const nlohmann::json doc =
        io::acceptanceRecordJson("acc_a", run, selection, flow, io::EventHalf::Odd, t);
    CHECK(doc.at("half") == "odd");
    CHECK(doc.at("source").at("run_id") == "run_a");
    CHECK_THAT(doc.at("entries").get<double>(), WithinULP(3.0, 0));

    const io::AcceptanceTable back =
        io::acceptanceTableFromRecord(nlohmann::json::parse(doc.dump()));
    REQUIRE(back.counts.size() == t.counts.size());
    for (std::size_t i = 0; i < t.counts.size(); ++i) {
        CHECK_THAT(back.counts[i], WithinULP(t.counts[i], 0));
    }
    CHECK_THAT(back.theta_edges.back(), WithinULP(kHalfPi, 0));

    nlohmann::json bad = doc;
    bad["acceptance_version"] = io::kAcceptanceRecordVersion + 1;
    CHECK_THROWS_AS(io::acceptanceTableFromRecord(bad), IoError);
    bad = doc;
    bad["counts"].erase(0);
    CHECK_THROWS_AS(io::acceptanceTableFromRecord(bad), IoError);
    bad = doc;
    bad.erase("En_edges_MeV");
    CHECK_THROWS_AS(io::acceptanceTableFromRecord(bad), IoError);
    CHECK_THAT(back.En_edges[1], WithinULP(11.0, 0));
}
