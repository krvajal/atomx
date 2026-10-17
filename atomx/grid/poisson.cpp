#include "atomx/grid/poisson.hpp"

#include <array>
#include <cmath>
#include <string>

namespace atomx::grid {

namespace {

// Solves the symmetric tridiagonal system with constant off-diagonal `off`
// and diagonal `diag`, overwriting `rhs` with the solution (Thomas algorithm).
void SolveTridiagonal(std::span<double> diag, double off,
                      std::span<double> rhs) {
  const std::size_t n = diag.size();
  for (std::size_t i = 1; i < n; ++i) {
    const double w = off / diag[i - 1];
    diag[i] -= w * off;
    rhs[i] -= w * rhs[i - 1];
  }
  rhs[n - 1] /= diag[n - 1];
  for (std::size_t i = n - 1; i-- > 0;) {
    rhs[i] = (rhs[i] - off * rhs[i + 1]) / diag[i];
  }
}

}  // namespace

Field RadialPoisson::Solve(int k, int nst, std::span<const double> f) const {
  const std::size_t mesh = grid_.size();
  const auto r = grid_.r();
  const auto r2 = grid_.r2();
  const auto sqr = grid_.sqrt_r();
  const double dx = grid_.dx();

  // Series of the solution near the origin.
  const int k21 = 2 * k + 1;
  const int nk1 = nst - k - 1;
  if (nk1 <= 0) {
    throw Error("RadialPoisson: need nst > k+1 (k=" + std::to_string(k) +
                ", nst=" + std::to_string(nst) + ")");
  }
  double c2 = 0.0, c3 = 0.0;
  if (nk1 < 3) {
    std::array<double, 4> d{};
    for (int i = 0; i < 4; ++i) d[i] = -k21 * f[i] / __builtin_powi(r[i], nst);
    const auto b = SmallRSeries(d, r);
    // The series of the solution starts at power nk1 of the source series.
    std::array<double, 5> e{};
    for (int i = 0; i < 4; ++i) e[nk1 - 1 + i] = b[i];
    c2 = e[0] / (4.0 * k + 6.0);
    c3 = e[1] / (6.0 * k + 12.0);
  }

  // Numerov discretisation in x = log r of (sqrt(r) * V).
  const double ch = dx * dx / 12.0;
  const double xkh2 = ch * (k + 0.5) * (k + 0.5);
  const double ei = 1.0 - xkh2;
  const double di = -(2.0 + 10.0 * xkh2);

  Field vh(mesh, 0.0);
  for (std::size_t i = 1; i < mesh; ++i) vh[i] = k21 * ch * sqr[i] * f[i];

  // Unknowns are points 1..mesh-2; the end points follow from the boundary
  // conditions and are eliminated from the first and last equations.
  const std::size_t m = mesh - 2;
  Field diag(m, -di);
  const double f1 = __builtin_powi(sqr[0] / sqr[1], k21);
  diag[0] -= ei * f1;
  const double fn = __builtin_powi(sqr[mesh - 2] / sqr[mesh - 1], k21);
  diag[m - 1] -= ei * fn;

  double vhim1 = k21 * ch * sqr[0] * f[0];
  for (std::size_t i = 1; i + 1 < mesh; ++i) {
    const double vhi = vh[i];
    vh[i] = vhim1 + 10.0 * vhi + vh[i + 1];
    vhim1 = vhi;
  }
  vh[1] -= ei * __builtin_powi(sqr[0], k21) *
           (c2 * (r2[1] - r2[0]) + c3 * (r[1] * r2[1] - r[0] * r2[0]));

  SolveTridiagonal(diag, -ei, std::span<double>(vh).subspan(1, m));

  const double c0 = vh[1] / __builtin_powi(sqr[1], k21) - c2 * r2[1] -
                    c3 * r[1] * r2[1];
  vh[0] = __builtin_powi(sqr[0], k21) * (c0 + c2 * r2[0] + c3 * r[0] * r2[0]);
  vh[mesh - 1] = vh[mesh - 2] * fn;

  for (std::size_t i = 0; i < mesh; ++i) vh[i] /= sqr[i];
  return vh;
}

}  // namespace atomx::grid
