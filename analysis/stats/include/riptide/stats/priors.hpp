#pragma once

// Proper priors of the Case A posterior (D33), returned as logarithms on a
// grid. They are plain functions passed to the model as arguments: the choice
// of prior is data, never a branch inside the model.
//
// Conventions: the E_n prior is a density in 1/MeV (it integrates to 1 over
// the domain); the theta_p priors are probability masses of the grid cells
// (they sum to 1), because theta_p is marginalised by a sum over the grid.

#include <vector>

#include "riptide/stats/grid.hpp"
#include "riptide/stats/numerics.hpp"

namespace riptide::stats {

// log(1 / (En_max - En_min)) at every point.
[[nodiscard]] std::vector<double> flatEnergyLogPrior(const Grid& En_grid);

// Isotropic CM scattering makes Ep uniform in [0, Ep_max]; in the
// non-relativistic limit Ep = En cos^2(theta_p), so theta_p has density
// sin(2 theta_p) on [0, pi/2]. The mass of each cell is the exact integral of
// that density over the cell, 0.5 (cos 2 lo - cos 2 hi), with the cell
// [theta - step/2, theta + step/2] clipped to [0, pi/2]: the end cells are
// half cells and the masses sum to 1 exactly up to rounding.
[[nodiscard]] std::vector<double> physicalRecoilAngleLogPrior(const Grid& theta_grid);

// -log(n) at every point: the prior of the linearised reference case R1.
[[nodiscard]] std::vector<double> flatRecoilAngleLogPrior(const Grid& theta_grid);

} // namespace riptide::stats
