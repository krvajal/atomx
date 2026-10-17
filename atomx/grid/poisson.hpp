#pragma once

#include <span>

#include "atomx/core/types.hpp"
#include "atomx/grid/radial_grid.hpp"

namespace atomx::grid {

// Radial multipole Poisson solver on a logarithmic grid (ld1.x's hartree):
//
//   V_k[f](r) = \int_0^\infty r_<^k / r_>^(k+1) f(r') dr'
//
// where f is a radial density (e.g. 4 pi r^2 n(r)) that behaves as r^nst near
// the origin. No e^2 factor is included. Solved as a Numerov tridiagonal
// system with series boundary conditions at the origin.
class RadialPoisson {
 public:
  explicit RadialPoisson(const RadialGrid& grid) : grid_(grid) {}

  Field Solve(int k, int nst, std::span<const double> f) const;

 private:
  const RadialGrid& grid_;
};

}  // namespace atomx::grid
