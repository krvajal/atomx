#pragma once

#include <memory>
#include <span>

#include "atomx/atom/nucleus.hpp"
#include "atomx/atom/orbital.hpp"
#include "atomx/core/types.hpp"
#include "atomx/grid/radial_grid.hpp"

namespace atomx::solvers {

struct RadialSolution {
  double energy = 0.0;  // Ry
  Spinor spinor;        // r * R(r), normalised
  bool converged = false;
};

// Bound-state eigensolver for one orbital in a local, spherical potential
// (Ry, including the nucleus). `energy_guess` is the previous eigenvalue.
class RadialSolver {
 public:
  virtual ~RadialSolver() = default;
  virtual RadialSolution Solve(const atom::Orbital& orbital,
                               std::span<const double> potential,
                               double energy_guess) const = 0;
};

// Schroedinger equation by Numerov integration with node counting and
// Cooley's energy correction (ld1.x's ascheq).
std::unique_ptr<RadialSolver> MakeSchrodingerSolver(
    const grid::RadialGrid& grid, const atom::Nucleus& nucleus);

// Dirac equation by predictor-corrector integration (ld1.x's dirsol).
// `speed_of_light` in Hartree atomic units.
std::unique_ptr<RadialSolver> MakeDiracSolver(const grid::RadialGrid& grid,
                                              const atom::Nucleus& nucleus,
                                              double speed_of_light);

}  // namespace atomx::solvers
