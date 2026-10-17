#pragma once

#include <span>

#include "atomx/core/types.hpp"
#include "atomx/grid/radial_grid.hpp"

namespace atomx::solvers {

// First-order response of a non-relativistic orbital (ld1.x's green):
// solves (H - e) y = -dvy for y orthogonal to chi, where H is the radial
// Schroedinger operator with potential v (Ry) and angular momentum l, chi
// the unperturbed orbital and e its eigenvalue. dvy must be orthogonal to chi.
// `ze2` = -Z e^2 of the point nucleus (series start at the origin).
class OrbitalResponse {
 public:
  OrbitalResponse(const grid::RadialGrid& grid, double ze2)
      : grid_(grid), ze2_(ze2) {}

  Field Solve(int l, double e, std::span<const double> dvy,
              std::span<const double> chi, std::span<const double> v) const;

 private:
  const grid::RadialGrid& grid_;
  double ze2_;
};

}  // namespace atomx::solvers
