#pragma once

#include <array>
#include <cstddef>
#include <span>

#include "atomx/core/types.hpp"

namespace atomx::grid {

// Coefficients b0..b3 of the cubic through the first four points of f:
// f(r) ~ b0 + b1 r + b2 r^2 + b3 r^3 near the origin.
std::array<double, 4> SmallRSeries(std::span<const double> f,
                                   std::span<const double> r);

// Logarithmic radial grid r_i = exp(xmin + i dx) / z, i = 0..n-1, with an odd
// number of points so that Simpson's rule closes. Same mesh as ld1.x's do_mesh.
class RadialGrid {
 public:
  static RadialGrid Logarithmic(double xmin, double dx, double rmax, double z);

  std::size_t size() const { return r_.size(); }
  double dx() const { return dx_; }
  double xmin() const { return xmin_; }

  double r(std::size_t i) const { return r_[i]; }
  std::span<const double> r() const { return r_; }
  std::span<const double> r2() const { return r2_; }
  std::span<const double> sqrt_r() const { return sqr_; }
  // dr/dx = r dx
  std::span<const double> rab() const { return rab_; }

  // \int_0^\infty f(r) dr over the first `n` points (default: all), with f
  // assumed to behave as r^nst near the origin (ld1.x's int_0_inf_dr).
  double Integrate(std::span<const double> f, int nst) const;
  double Integrate(std::span<const double> f, int nst, std::size_t n) const;

  Field Zeros() const { return Field(size(), 0.0); }

 private:
  RadialGrid() = default;

  double xmin_ = 0.0;
  double dx_ = 0.0;
  Field r_, r2_, sqr_, rab_;
};

}  // namespace atomx::grid
