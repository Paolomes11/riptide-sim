#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/stats/acceptance.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"

namespace st = riptide::stats;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace {

constexpr double kHalfPi = 0.5 * riptide::units::pi;
constexpr std::size_t kThetaBins = 18;
constexpr std::size_t kEnBins = 10;
constexpr double kEnLo = 10.0;
constexpr double kEnHi = 19.5;
constexpr double kEnStep = 1.0;
constexpr double kEvents = 1.0e5;

std::vector<double> linspace(double lo, double hi, std::size_t bins) {
    std::vector<double> edges(bins + 1);
    for (std::size_t i = 0; i <= bins; ++i) {
        edges[i] = lo + ((hi - lo) * static_cast<double>(i) / static_cast<double>(bins));
    }
    edges.back() = hi;
    return edges;
}

// 1 MeV E_n bins from 10 MeV with a last partial bin up to 19.5 MeV, the
// binning of D40.
std::vector<double> energyEdges() {
    std::vector<double> edges;
    for (std::size_t k = 0; k < kEnBins; ++k) {
        edges.push_back(kEnLo + (kEnStep * static_cast<double>(k)));
    }
    edges.push_back(kEnHi);
    return edges;
}

// Counts n_kt = shape(k, t) * m_iso(t; E_k) * kEvents: the ratio of bin
// (k, t) is then shape(k, t) / sum_t shape(k, t) m_iso(t; E_k).
template <typename Shape> st::AcceptanceHistogram histogramWithShape(Shape shape) {
    st::AcceptanceHistogram h;
    h.En_edges = energyEdges();
    h.theta_edges = linspace(0.0, kHalfPi, kThetaBins);
    for (std::size_t k = 0; k < kEnBins; ++k) {
        const double En = 0.5 * (h.En_edges[k] + h.En_edges[k + 1]);
        for (std::size_t t = 0; t < kThetaBins; ++t) {
            const double iso =
                st::isotropicRecoilAngleMass(En, h.theta_edges[t], h.theta_edges[t + 1]);
            h.counts.push_back(shape(k, t) * iso * kEvents);
        }
    }
    return h;
}

double one(std::size_t /*k*/, std::size_t /*t*/) {
    return 1.0;
}

} // namespace

TEST_CASE("the isotropic theta_p mass is relativistic and sums to 1", "[stats][acceptance]") {
    for (const double En : {1.0, 10.0, 14.0, 20.0, 30.0}) {
        CHECK_THAT(st::isotropicRecoilAngleMass(En, 0.0, kHalfPi), WithinAbs(1.0, 1e-15));
        const std::vector<double> edges = linspace(0.0, kHalfPi, 37);
        double sum = 0.0;
        for (std::size_t t = 0; t + 1 < edges.size(); ++t) {
            const double m = st::isotropicRecoilAngleMass(En, edges[t], edges[t + 1]);
            CHECK(m > 0.0);
            // Close to the non-relativistic 0.5 (cos 2lo - cos 2hi): the
            // difference is ~1% at 20 MeV.
            const double nonrel = 0.5 * (std::cos(2.0 * edges[t]) - std::cos(2.0 * edges[t + 1]));
            CHECK_THAT(m, WithinRel(nonrel, 0.02));
            sum += m;
        }
        CHECK_THAT(sum, WithinAbs(1.0, 1e-14));
    }
    CHECK_THAT(st::isotropicRecoilAngleMass(14.0, 0.3, 0.3), WithinAbs(0.0, 0.0));
    CHECK_THROWS_AS(st::isotropicRecoilAngleMass(0.0, 0.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(st::isotropicRecoilAngleMass(14.0, 0.5, 0.4), std::invalid_argument);
    CHECK_THROWS_AS(st::isotropicRecoilAngleMass(14.0, -0.1, 0.4), std::invalid_argument);
    CHECK_THROWS_AS(st::isotropicRecoilAngleMass(14.0, 0.0, 1.6), std::invalid_argument);
}

TEST_CASE("isotropic counts give r = 1 and the isotropic prior", "[stats][acceptance]") {
    const st::AcceptanceRatio ratio(histogramWithShape(one));
    for (const double En : {5.0, 10.0, 14.3, 19.25, 25.0}) {
        for (const double theta : {0.0, 0.01, 0.7, 1.2, kHalfPi}) {
            CHECK_THAT(ratio(En, theta), WithinAbs(1.0, 1e-12));
        }
    }
    CHECK_THAT(ratio.maximum(), WithinAbs(1.0, 1e-12));

    const st::Grid En = st::energyGrid(st::EnergyDomain{});
    const st::Grid theta = st::recoilAngleGrid();
    const auto rows = st::acceptanceRecoilAngleLogPrior(ratio, En, theta);
    REQUIRE(rows.size() == En.size());
    const double half = 0.5 * theta.step();
    for (std::size_t i = 0; i < En.size(); i += 37) {
        REQUIRE(rows[i].size() == theta.size());
        for (std::size_t j = 0; j < theta.size(); ++j) {
            const double lo = std::max(0.0, theta[j] - half);
            const double hi = std::min(kHalfPi, theta[j] + half);
            const double expected = std::log(st::isotropicRecoilAngleMass(En[i], lo, hi));
            CHECK_THAT(rows[i][j], WithinAbs(expected, 1e-11));
        }
    }
}

TEST_CASE("the acceptance ratio interpolates bilinearly and is flat outside",
          "[stats][acceptance]") {
    // shape = 1 + k + t / 2: linear in both bin indices, so linear between
    // the centres up to the per-row normalisation of the ratio.
    const auto shape = [](std::size_t k, std::size_t t) {
        return 1.0 + static_cast<double>(k) + (0.5 * static_cast<double>(t));
    };
    const st::AcceptanceRatio ratio(histogramWithShape(shape));
    const auto& Ec = ratio.energyCentres();
    const auto& tc = ratio.angleCentres();
    REQUIRE(Ec.size() == kEnBins);
    REQUIRE(tc.size() == kThetaBins);
    CHECK_THAT(Ec.back(), WithinULP(19.25, 0));

    // Within one E_n row, nodes are proportional to the shape.
    const double r00 = ratio(Ec[0], tc[0]);
    CHECK_THAT(ratio(Ec[0], tc[4]) / r00, WithinRel(shape(0, 4) / shape(0, 0), 1e-12));
    // Linear between two theta nodes.
    const double mid = 0.5 * (tc[2] + tc[3]);
    CHECK_THAT(ratio(Ec[5], mid),
               WithinRel(0.5 * (ratio(Ec[5], tc[2]) + ratio(Ec[5], tc[3])), 1e-12));
    // Linear between two E_n nodes.
    const double Emid = 0.5 * (Ec[3] + Ec[4]);
    CHECK_THAT(ratio(Emid, tc[7]),
               WithinRel(0.5 * (ratio(Ec[3], tc[7]) + ratio(Ec[4], tc[7])), 1e-12));
    // Constant beyond the outer centres, in both variables.
    CHECK_THAT(ratio(9.0, tc[6]), WithinULP(ratio(Ec.front(), tc[6]), 0));
    CHECK_THAT(ratio(30.0, tc[6]), WithinULP(ratio(Ec.back(), tc[6]), 0));
    CHECK_THAT(ratio(Ec[2], 0.0), WithinULP(ratio(Ec[2], tc.front()), 0));
    CHECK_THAT(ratio(Ec[2], kHalfPi), WithinULP(ratio(Ec[2], tc.back()), 0));
    // The maximum is the largest node.
    double largest = 0.0;
    for (const double E : Ec) {
        for (const double t : tc) {
            largest = std::max(largest, ratio(E, t));
        }
    }
    CHECK_THAT(ratio.maximum(), WithinULP(largest, 0));
}

TEST_CASE("acceptance prior rows are normalised and follow r", "[stats][acceptance]") {
    // Forward deficit like the measured one: half the isotropic rate in the
    // first bin, rising to 1 by 30 deg.
    const auto shape = [](std::size_t /*k*/, std::size_t t) {
        return t < 6 ? 0.5 + (0.5 * static_cast<double>(t) / 6.0) : 1.0;
    };
    const st::AcceptanceRatio ratio(histogramWithShape(shape));
    const st::Grid En = st::energyGrid(st::EnergyDomain{0.5 * kEnLo, 3.0 * kEnLo});
    const st::Grid theta = st::recoilAngleGrid();
    const auto rows = st::acceptanceRecoilAngleLogPrior(ratio, En, theta);
    const std::vector<double> iso = st::physicalRecoilAngleLogPrior(theta);
    for (const auto& row : rows) {
        double sum = 0.0;
        for (const double lp : row) {
            REQUIRE(std::isfinite(lp));
            sum += std::exp(lp);
        }
        CHECK_THAT(sum, WithinAbs(1.0, 1e-12));
        // Forward cells lose mass relative to the isotropic prior, the
        // large-angle ones gain it.
        CHECK(row.front() < iso.front());
        CHECK(row.back() > iso.back());
    }
}

TEST_CASE("invalid acceptance histograms are rejected", "[stats][acceptance]") {
    const st::AcceptanceHistogram good = histogramWithShape(one);
    CHECK_NOTHROW(st::AcceptanceRatio(good));

    auto bad = good;
    bad.counts[7] = 0.0;
    CHECK_THROWS_AS(st::AcceptanceRatio(bad), std::invalid_argument);
    bad = good;
    bad.counts[3] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(st::AcceptanceRatio(bad), std::invalid_argument);
    bad = good;
    bad.counts.pop_back();
    CHECK_THROWS_AS(st::AcceptanceRatio(bad), std::invalid_argument);
    bad = good;
    bad.En_edges[3] = bad.En_edges[2];
    CHECK_THROWS_AS(st::AcceptanceRatio(bad), std::invalid_argument);
    bad = good;
    bad.En_edges.front() = 0.0;
    CHECK_THROWS_AS(st::AcceptanceRatio(bad), std::invalid_argument);
    bad = good;
    bad.theta_edges.back() = 1.6;
    CHECK_THROWS_AS(st::AcceptanceRatio(bad), std::invalid_argument);
    bad = good;
    bad.theta_edges = {0.0};
    bad.counts.clear();
    CHECK_THROWS_AS(st::AcceptanceRatio(bad), std::invalid_argument);
}
