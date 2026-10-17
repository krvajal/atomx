#include "atomx/grid/radial_grid.hpp"

#include <cmath>

namespace atomx::grid {

std::array<double, 4> SmallRSeries(std::span<const double> f,
                                   std::span<const double> r) {
  const double dr21 = r[1] - r[0];
  const double dr31 = r[2] - r[0];
  const double dr32 = r[2] - r[1];
  const double dr41 = r[3] - r[0];
  const double dr42 = r[3] - r[1];
  const double dr43 = r[3] - r[2];
  const double df21 = (f[1] - f[0]) / dr21;
  const double df32 = (f[2] - f[1]) / dr32;
  const double df43 = (f[3] - f[2]) / dr43;
  const double ddf42 = (df43 - df32) / dr42;
  const double ddf31 = (df32 - df21) / dr31;
  std::array<double, 4> b{};
  b[3] = (ddf42 - ddf31) / dr41;
  b[2] = ddf31 - b[3] * (r[0] + r[1] + r[2]);
  b[1] = df21 - b[2] * (r[1] + r[0]) -
         b[3] * (r[0] * r[0] + r[1] * r[1] + r[0] * r[1]);
  b[0] = f[0] - r[0] * (b[1] + r[0] * (b[2] + r[0] * b[3]));
  return b;
}

RadialGrid RadialGrid::Logarithmic(double xmin, double dx, double rmax,
                                   double z) {
  if (dx <= 0.0 || z <= 0.0 || rmax <= 0.0) {
    throw Error("RadialGrid: dx, z and rmax must be positive");
  }
  const double xmax = std::log(rmax * z);
  auto mesh = static_cast<std::size_t>((xmax - xmin) / dx + 1.0);
  mesh = 2 * (mesh / 2) + 1;  // odd, for Simpson's rule

  RadialGrid g;
  g.xmin_ = xmin;
  g.dx_ = dx;
  g.r_.resize(mesh);
  g.r2_.resize(mesh);
  g.sqr_.resize(mesh);
  g.rab_.resize(mesh);
  for (std::size_t i = 0; i < mesh; ++i) {
    const double x = xmin + static_cast<double>(i) * dx;
    g.r_[i] = std::exp(x) / z;
    g.r2_[i] = g.r_[i] * g.r_[i];
    g.sqr_[i] = std::sqrt(g.r_[i]);
    g.rab_[i] = g.r_[i] * dx;
  }
  return g;
}

double RadialGrid::Integrate(std::span<const double> f, int nst) const {
  return Integrate(f, nst, size());
}

double RadialGrid::Integrate(std::span<const double> f, int nst,
                             std::size_t n) const {
  // Analytic piece from 0 to r_0, using the series of f / r^nst.
  std::array<double, 4> fs{};
  for (int i = 0; i < 4; ++i) fs[i] = f[i] / __builtin_powi(r_[i], nst);
  const auto b = SmallRSeries(fs, r_);
  const double r0 = r_[0];
  double sum = (b[0] / (nst + 1) +
                r0 * (b[1] / (nst + 2) + r0 * b[2] / (nst + 3))) *
               __builtin_powi(r0, nst + 1);
  // Simpson's rule in x, dr = r dx.
  double s = 0.0;
  for (std::size_t i = 0; i + 2 < n; i += 2) {
    s += f[i] * r_[i] + 4.0 * f[i + 1] * r_[i + 1] + f[i + 2] * r_[i + 2];
  }
  return sum + s * dx_ / 3.0;
}

}  // namespace atomx::grid
