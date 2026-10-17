// Port of ld1.x's ascheq + start_scheq (Quantum ESPRESSO, GPL).
// Internal arrays are 1-based (index 0 unused) to follow the original.

#include <cmath>
#include <vector>

#include "atomx/core/constants.hpp"
#include "atomx/grid/radial_grid.hpp"
#include "atomx/solvers/radial_solver.hpp"

namespace atomx::solvers {

namespace {

class SchrodingerSolver final : public RadialSolver {
 public:
  SchrodingerSolver(const grid::RadialGrid& grid, const atom::Nucleus& nucleus)
      : grid_(grid), ze2_(-nucleus.charge() * constants::kE2) {}

  RadialSolution Solve(const atom::Orbital& orbital,
                       std::span<const double> potential,
                       double energy_guess) const override;

 private:
  static constexpr int kMaxIterations = 50;
  static constexpr double kThreshold = 1.0e-10;

  const grid::RadialGrid& grid_;
  double ze2_;
};

double SignOf(double a, double b) { return b >= 0.0 ? std::abs(a) : -std::abs(a); }

RadialSolution SchrodingerSolver::Solve(const atom::Orbital& orbital,
                                        std::span<const double> vpot0,
                                        double e) const {
  const int mesh = static_cast<int>(grid_.size());
  const int lam = orbital.l;
  const int nn = orbital.n;
  // 1-based views
  std::vector<double> r(mesh + 1), r2(mesh + 1), sqr(mesh + 1), v(mesh + 1);
  for (int i = 1; i <= mesh; ++i) {
    r[i] = grid_.r(i - 1);
    r2[i] = grid_.r2()[i - 1];
    sqr[i] = grid_.sqrt_r()[i - 1];
    v[i] = vpot0[i - 1];
  }
  std::vector<double> y(mesh + 1, 0.0), f(mesh + 1), c(mesh + 1),
      el(mesh + 1);
  const double dx = grid_.dx();
  const double ze2 = ze2_;

  RadialSolution out;
  // ld1.x keeps the last trial energy when an orbital fails to converge.
  auto Failed = [&](double energy) {
    out.energy = energy;
    out.spinor = Spinor{Field(mesh, 0.0), Field(mesh, 0.0)};
    return out;
  };
  double thresh = kThreshold;
  if (e < -5.0e2) thresh *= 10.0;
  const double ddx12 = dx * dx / 12.0;
  const int l1 = lam + 1;
  const double sqlhf = (lam + 0.5) * (lam + 0.5);
  const int ndcr = nn - lam - 1;

  double eup = v[mesh] + sqlhf / r2[mesh];
  double elw = eup;
  for (int i = 1; i <= mesh; ++i) elw = std::min(elw, v[i] + sqlhf / r2[i]);
  if (eup == elw) return Failed(e);
  if (e > eup) e = 0.9 * eup + 0.1 * elw;
  if (e < elw) e = 0.9 * elw + 0.1 * eup;

  // Series of the potential near the origin, without the nuclear term.
  std::array<double, 4> ys{};
  for (int i = 1; i <= 4; ++i) ys[i - 1] = v[i] - ze2 / r[i];
  const auto b = grid::SmallRSeries(ys, std::span<const double>(&r[1], 4));

  bool converged = false;
  int ik = 0, nstart = mesh;
  double sum = 0.0, f2 = 0.0;
  for (int iter = 1;; ++iter) {
    if (iter > kMaxIterations) return Failed(e);

    f[1] = ddx12 * (r2[1] * (v[1] - e) + sqlhf);
    for (int i = 2; i <= mesh; ++i) {
      f[i] = ddx12 * (r2[i] * (v[i] - e) + sqlhf);
      if (f[i] != SignOf(f[i], f[i - 1])) ik = i;
    }
    if (ik >= mesh - 2) return Failed(e);
    for (int i = 1; i <= mesh; ++i) f[i] = 1.0 - f[i];
    std::fill(y.begin(), y.end(), 0.0);

    // First two points from the series (start_scheq).
    const double xl1 = lam + 1.0;
    const double b0e = b[0] - e;
    const double c1 = 0.5 * ze2 / xl1;
    const double c2 = (c1 * ze2 + b0e) / (4.0 * lam + 6.0);
    const double c3 = (c2 * ze2 + c1 * b0e + b[1]) / (6.0 * lam + 12.0);
    const double c4 =
        (c3 * ze2 + c2 * b0e + c1 * b[1] + b[2]) / (8.0 * lam + 20.0);
    for (int i = 1; i <= 2; ++i) {
      const double x = r[i];
      y[i] = (1.0 + x * (c1 + x * (c2 + x * (c3 + x * c4)))) *
             __builtin_powi(x, l1) / sqr[i];
    }

    // Outward integration to the matching point, counting nodes.
    int ncross = 0;
    double ymx = 0.0;
    for (int n = 2; n <= ik - 1; ++n) {
      y[n + 1] = ((12.0 - 10.0 * f[n]) * y[n] - f[n - 1] * y[n - 1]) / f[n + 1];
      if (y[n] != SignOf(y[n], y[n + 1])) ++ncross;
      ymx = std::max(ymx, std::abs(y[n + 1]));
    }
    if (ndcr < ncross) {
      eup = e;
      const double rap = std::pow(double(ncross + l1) / nn, 2);
      e = (e - v[mesh]) * rap + v[mesh];
      if (e < elw) e = 0.9 * elw + 0.1 * eup;
      continue;
    }
    if (ndcr > ncross) {
      elw = e;
      const double rap = std::pow(double(ncross + l1) / nn, 2);
      e = (e - v[mesh]) * rap + v[mesh];
      if (e > eup) e = 0.9 * eup + 0.1 * elw;
      continue;
    }

    // Inward integration (Froese's tail procedure).
    nstart = mesh;
    const double rstart = 10.0 * r[ik];
    if (rstart < r[mesh]) {
      for (int i = ik; i <= mesh; ++i) {
        nstart = i;
        if (r[i] >= rstart) break;
      }
      nstart = 2 * (nstart / 2) + 1;
    }
    el[ik + 1] = 10.0 * f[ik + 1] - 12.0;
    c[ik + 1] = -f[ik] * y[ik];
    for (int n = ik + 2; n <= nstart; ++n) {
      const double di = 10.0 * f[n] - 12.0;
      el[n] = di - f[n] * f[n - 1] / el[n - 1];
      c[n] = -c[n - 1] * f[n - 1] / el[n - 1];
    }
    double expn = std::exp(-std::sqrt(12.0 * std::abs(1.0 - f[nstart - 1])));
    y[nstart - 1] = c[nstart - 1] / (el[nstart - 1] + f[nstart] * expn);
    y[nstart] = expn * y[nstart - 1];
    for (int n = nstart - 2; n >= ik + 1; --n) {
      y[n] = (c[n] - f[n + 1] * y[n + 1]) / el[n];
    }

    // Cooley's correction of the eigenvalue.
    const double fe =
        (12.0 - 10.0 * f[ik]) * y[ik] - f[ik - 1] * y[ik - 1] - f[ik + 1] * y[ik + 1];
    if (ymx >= 1.0e10) {
      for (int i = 1; i <= mesh; ++i) y[i] /= ymx;
    }
    const double a0 = 1.0 / (2 * lam + 3);
    const double a1 = c1 / (lam + 2);
    const double a2 = (c1 * c1 + c2 + c2) / (2 * lam + 5);
    const double sum0 =
        (a0 + r[1] * (a1 + r[1] * a2)) * __builtin_powi(r[1], 2 * lam + 3);
    f2 = r2[1] * y[1] * y[1];
    sum = r[1] * f2 / (2 * l1 + 1);
    for (int n = 1; n <= nstart - 2; n += 2) {
      const double f0 = f2;
      const double f1 = r2[n + 1] * y[n + 1] * y[n + 1];
      f2 = r2[n + 2] * y[n + 2] * y[n + 2];
      sum += f0 + f2 + 4.0 * f1;
    }
    sum = sum0 + dx * sum / 3.0;
    const double dfe = -y[ik] * f[ik] / dx / sum;
    double de = -fe * dfe;
    const double eps = std::abs(de / e);
    if (std::abs(de) < thresh) {
      converged = true;
      break;
    }
    if (eps > 0.25) de = 0.25 * de / eps;
    if (de > 0.0) elw = e;
    if (de < 0.0) eup = e;
    e += de;
    if (e > eup) e = 0.9 * eup + 0.1 * elw;
    if (e < elw) e = 0.9 * elw + 0.1 * eup;
    if (iter >= kMaxIterations) break;  // normalise anyway, flag failure
  }

  // Exponential tail and normalisation.
  for (int n = nstart; n <= mesh - 1; ++n) {
    y[n + 1] = 0.0;
    if (y[n] == 0.0) continue;
    const double expn = std::log(std::abs(y[n])) -
                        std::sqrt(12.0 * std::abs(1.0 - f[n]));
    if (expn < -80.0) continue;
    y[n + 1] = SignOf(std::exp(expn), y[n]);
  }
  double sum1 = 0.0;
  for (int n = nstart; n <= mesh - 2; n += 2) {
    const double f0 = f2;
    const double f1 = r2[n + 1] * y[n + 1] * y[n + 1];
    f2 = r2[n + 2] * y[n + 2] * y[n + 2];
    sum1 += f0 + f2 + 4.0 * f1;
  }
  sum = std::sqrt(sum + dx * sum1 / 3.0);

  out.energy = e;
  out.converged = converged;
  out.spinor.large.resize(mesh);
  out.spinor.small.assign(mesh, 0.0);
  for (int n = 1; n <= mesh; ++n) out.spinor.large[n - 1] = sqr[n] * y[n] / sum;
  return out;
}

}  // namespace

std::unique_ptr<RadialSolver> MakeSchrodingerSolver(
    const grid::RadialGrid& grid, const atom::Nucleus& nucleus) {
  return std::make_unique<SchrodingerSolver>(grid, nucleus);
}

}  // namespace atomx::solvers
