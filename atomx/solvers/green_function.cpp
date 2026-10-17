// Port of ld1.x's green + outward + inward (Quantum ESPRESSO, GPL).
// Internal arrays are 1-based (index 0 unused) to follow the original.

#include "atomx/solvers/green_function.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace atomx::solvers {

Field OrbitalResponse::Solve(int lam, double e, std::span<const double> dvy0,
                             std::span<const double> chi,
                             std::span<const double> v0) const {
  const int mesh = static_cast<int>(grid_.size());
  std::vector<double> r(mesh + 1), r2(mesh + 1), sqr(mesh + 1), v(mesh + 1),
      dvy(mesh + 1), f(mesh + 1), g(mesh + 1), y(mesh + 1, 0.0);
  for (int i = 1; i <= mesh; ++i) {
    r[i] = grid_.r(i - 1);
    r2[i] = grid_.r2()[i - 1];
    sqr[i] = grid_.sqrt_r()[i - 1];
    v[i] = v0[i - 1];
    dvy[i] = dvy0[i - 1];
  }
  const double ze2 = ze2_;
  const double ddx12 = grid_.dx() * grid_.dx() / 12.0;
  const int l1 = lam + 1;
  const double sqlhf = (lam + 0.5) * (lam + 0.5);
  const double xl1 = l1;

  // Series of the potential (without the nucleus) and of dvy / r^(l+1).
  std::array<double, 4> tmp{};
  for (int i = 1; i <= 4; ++i) tmp[i - 1] = v[i] - ze2 / r[i];
  const auto b = grid::SmallRSeries(tmp, std::span<const double>(&r[1], 4));
  for (int i = 1; i <= 4; ++i) tmp[i - 1] = dvy[i] / __builtin_powi(r[i], l1);
  const auto a = grid::SmallRSeries(tmp, std::span<const double>(&r[1], 4));

  int imatch = 1;
  f[1] = ddx12 * (r2[1] * (v[1] - e) + sqlhf);
  for (int i = 2; i <= mesh; ++i) {
    f[i] = ddx12 * (r2[i] * (v[i] - e) + sqlhf);
    // Fortran: f(i) /= sign(f(i), f(i-1))
    if (f[i] != (f[i - 1] >= 0.0 ? std::abs(f[i]) : -std::abs(f[i]))) imatch = i;
  }
  for (int i = 1; i <= mesh; ++i) {
    f[i] = 1.0 - f[i];
    g[i] = ddx12 * dvy[i] * r[i] * sqr[i];
  }

  const double b0e = b[0] - e;
  const double c0 = 1.0;
  const double c1 = 0.5 * ze2 * c0 / xl1;
  const double c2 = (c1 * ze2 + b0e * c0 + a[0]) / (4 * lam + 6);
  const double c3 = (c2 * ze2 + c1 * b0e + b[1] * c0 + a[1]) / (6 * lam + 12);
  const double c4 =
      (c3 * ze2 + c2 * b0e + c1 * b[1] + b[2] * c0 + a[2]) / (8 * lam + 20);
  for (int i = 1; i <= 2; ++i) {
    const double x = r[i];
    y[i] = (c0 + x * (c1 + x * (c2 + x * (c3 + x * c4)))) * __builtin_powi(x, l1) /
           sqr[i];
  }
  // Numerov form of the inhomogeneous term.
  double gim1 = g[1];
  for (int i = 2; i <= mesh - 1; ++i) {
    const double gi = g[i];
    g[i] = gim1 + 10.0 * gi + g[i + 1];
    gim1 = gi;
  }

  // Outward to the matching point.
  for (int n = 2; n <= imatch - 1; ++n) {
    y[n + 1] =
        ((12.0 - 10.0 * f[n]) * y[n] - f[n - 1] * y[n - 1] + g[n]) / f[n + 1];
  }
  // Inward from the end of the mesh.
  {
    std::vector<double> el(mesh + 1), c(mesh + 1);
    const int imp1 = imatch + 1, nm1 = mesh - 1;
    el[imp1] = 10.0 * f[imp1] - 12.0;
    c[imp1] = g[imp1] - f[imatch] * y[imatch];
    for (int n = imatch + 2; n <= nm1; ++n) {
      const double di = 10.0 * f[n] - 12.0;
      el[n] = di - f[n] * f[n - 1] / el[n - 1];
      c[n] = g[n] - c[n - 1] * f[n - 1] / el[n - 1];
    }
    const double ff = 12.0 * std::abs(1.0 - f[nm1]);
    const double expn = std::exp(-std::sqrt(ff));
    y[mesh] = (c[nm1] * expn + g[nm1] * el[nm1] * (expn - 1.0) / ff) /
              (el[nm1] + f[mesh] * expn);
    for (int n = nm1; n >= imp1; --n) y[n] = (c[n] - f[n + 1] * y[n + 1]) / el[n];
  }

  Field out(mesh), work(mesh);
  for (int i = 1; i <= mesh; ++i) {
    out[i - 1] = y[i] * sqr[i];
    work[i - 1] = out[i - 1] * chi[i - 1];
  }
  const double fac = grid_.Integrate(work, 2 * lam + 2);
  for (int i = 0; i < mesh; ++i) out[i] -= fac * chi[i];
  return out;
}

}  // namespace atomx::solvers
