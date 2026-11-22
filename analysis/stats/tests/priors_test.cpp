#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"

namespace st = riptide::stats;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

double trapezoid(const st::Grid& grid, const std::vector<double>& log_density) {
    double sum = 0.0;
    for (std::size_t i = 1; i < grid.size(); ++i) {
        sum += 0.5 * (std::exp(log_density[i - 1]) + std::exp(log_density[i])) * grid.step();
    }
    return sum;
}

double massSum(const std::vector<double>& log_mass) {
    double sum = 0.0;
    for (const double value : log_mass) {
        sum += std::exp(value);
    }
    return sum;
}

} // namespace

TEST_CASE("the flat E_n prior is a proper density", "[stats][priors]") {
    const st::Grid En = st::energyGrid(st::EnergyDomain{});
    const std::vector<double> log_prior = st::flatEnergyLogPrior(En);
    REQUIRE(log_prior.size() == En.size());
    CHECK_THAT(trapezoid(En, log_prior), WithinRel(1.0, 1e-12));
    CHECK_THAT(std::exp(log_prior.front()), WithinRel(0.1, 1e-12)); // 1 / (10 MeV)
}

TEST_CASE("the sin(2 theta) prior has exact cell masses", "[stats][priors]") {
    const st::Grid theta = st::recoilAngleGrid();
    const std::vector<double> log_mass = st::physicalRecoilAngleLogPrior(theta);
    REQUIRE(log_mass.size() == theta.size());
    CHECK_THAT(massSum(log_mass), WithinAbs(1.0, 1e-12));
    // Inside, the exact mass is sin(2 theta) * step up to O(step^2): the same
    // tolerances as the toy test.
    for (std::size_t j = 1; j + 1 < theta.size(); ++j) {
        const double approx = std::sin(2.0 * theta[j]) * theta.step();
        CHECK_THAT(std::exp(log_mass[j]), WithinRel(approx, 1e-4) || WithinAbs(approx, 1e-9));
    }
    // The end cells are half cells where the density vanishes, so they are
    // the two smallest masses but still positive: no -inf in the prior.
    CHECK(std::isfinite(log_mass.front()));
    CHECK(std::isfinite(log_mass.back()));
    CHECK(log_mass.front() < log_mass[1]);
    CHECK(log_mass.back() < log_mass[theta.size() - 2]);
}

TEST_CASE("the flat theta prior is uniform over the cells", "[stats][priors]") {
    const st::Grid theta = st::recoilAngleGrid();
    const std::vector<double> log_mass = st::flatRecoilAngleLogPrior(theta);
    CHECK_THAT(massSum(log_mass), WithinAbs(1.0, 1e-12));
    CHECK_THAT(log_mass[7], WithinRel(-std::log(500.0), 1e-15));
}
