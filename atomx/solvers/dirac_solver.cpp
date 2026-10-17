// Port of ld1.x's dirsol + cfdsol + nodeno (Quantum ESPRESSO, GPL).
//
// Radial Dirac equation in Rydberg units, f = r * minor, g = r * major:
//   df/dr = +k/r f - (e - v) alpha/2 g
//   dg/dr = -k/r g + (e - v + 4/alpha^2) alpha/2 f
// Internal arrays are 1-based (index 0 unused) to follow the original.

#include <array>
#include <cmath>
#include <vector>

#include "atomx/core/constants.hpp"
#include "atomx/grid/radial_grid.hpp"
#include "atomx/solvers/radial_solver.hpp"

namespace atomx::solvers {

namespace {

// 2x2 coefficient matrix of the ODE at every grid point, in x = log r.
struct Coefficients {
  std::vector<std::array<std::array<double, 2>, 2>> z;  // z[i][a][b]
};

// Fifth-order Adams predictor-corrector from point j1 to j2 (either
// direction); y must hold 5 starting values before j1 (ld1.x's cfdsol).
void PredictorCorrector(const Coefficients& zz,
                        std::vector<std::array<double, 2>>& yy, int j1,
                        int j2) {
  static constexpr double abp[6] = {0.0,           1901.0 / 720.0,
                                    -1387.0 / 360.0, 109.0 / 30.0,
                                    -637.0 / 360.0,  251.0 / 720.0};
  static constexpr double amc[5] = {251.0 / 720.0, 323.0 / 360.0,
                                    -11.0 / 30.0, 53.0 / 360.0,
                                    -19.0 / 720.0};
  const int isgn = (j2 - j1) > 0 ? 1 : -1;
  const auto& z = zz.z;
  double fa[6], fb[6];
  for (int j = 1; j <= 5; ++j) {
    const int ip = j1 - isgn * j;
    fa[j] = z[ip][0][0] * yy[ip][0] + z[ip][0][1] * yy[ip][1];
    fb[j] = z[ip][1][0] * yy[ip][0] + z[ip][1][1] * yy[ip][1];
  }
  for (int j = j1; j != j2 + isgn; j += isgn) {
    double arp = yy[j - isgn][0];
    double brp = yy[j - isgn][1];
    for (int i = 1; i <= 5; ++i) {
      arp += isgn * abp[i] * fa[i];
      brp += isgn * abp[i] * fb[i];
    }
    fa[0] = z[j][0][0] * arp + z[j][0][1] * brp;
    fb[0] = z[j][1][0] * arp + z[j][1][1] * brp;
    yy[j][0] = yy[j - isgn][0];
    yy[j][1] = yy[j - isgn][1];
    for (int i = 0; i <= 4; ++i) {
      yy[j][0] += isgn * amc[i] * fa[i];
      yy[j][1] += isgn * amc[i] * fb[i];
    }
    for (int i = 5; i >= 2; --i) {
      fa[i] = fa[i - 1];
      fb[i] = fb[i - 1];
    }
    fa[1] = z[j][0][0] * yy[j][0] + z[j][0][1] * yy[j][1];
    fb[1] = z[j][1][0] * yy[j][0] + z[j][1][1] * yy[j][1];
  }
}

class DiracSolver final : public RadialSolver {
 public:
  DiracSolver(const grid::RadialGrid& grid, const atom::Nucleus& nucleus,
              double speed_of_light)
      : grid_(grid),
        z_(nucleus.charge()),
        point_nucleus_(nucleus.IsPointLike()),
        tbya_(2.0 * speed_of_light),
        abyt_(1.0 / (2.0 * speed_of_light)) {}

  RadialSolution Solve(const atom::Orbital& orbital,
                       std::span<const double> potential,
                       double energy_guess) const override;

 private:
  static constexpr int kMaxIterations = 100;
  static constexpr double kThreshold = 1.0e-10;

  // Values of (f, g) at the first five points from the series at r -> 0.
  void StartAtOrigin(int kcur, double ecur, double v1,
                     std::vector<std::array<double, 2>>& yy, int& ig) const;

  const grid::RadialGrid& grid_;
  double z_;
  bool point_nucleus_;
  double tbya_;  // 2 / alpha
  double abyt_;  // alpha / 2
};

void DiracSolver::StartAtOrigin(int kcur, double ecur, double v1,
                                std::vector<std::array<double, 2>>& yy,
                                int& ig) const {
  const auto r = grid_.r();
  const double tbya = tbya_, abyt = abyt_;
  double f0, g0, power;
  double vzero;
  if (point_nucleus_) {
    // Coulomb singularity: f, g ~ r^gamma, gamma = sqrt(k^2 - (alpha Z)^2).
    vzero = v1 + z_ * 2.0 / r[0];
    const double gamma0 = std::sqrt(kcur * kcur - 4.0 * (abyt * z_) * (abyt * z_));
    if (kcur < 0) {
      ig = -kcur;
      f0 = (kcur + gamma0) / (2.0 * abyt * z_);
      g0 = 1.0;
    } else {
      ig = kcur;
      f0 = 1.0;
      g0 = (kcur - gamma0) / (2.0 * abyt * z_);
    }
    power = gamma0;
    const double f1 = -(ecur - vzero) * abyt * g0 / (gamma0 - kcur + 1.0);
    const double g1 = (ecur - vzero + tbya * tbya) * abyt * f0 / (gamma0 + kcur + 1.0);
    const double f2 = -(ecur - vzero) * abyt * g1 / (gamma0 - kcur + 2.0);
    const double g2 = (ecur - vzero + tbya * tbya) * abyt * f1 / (gamma0 + kcur + 2.0);
    for (int ir = 1; ir <= 5; ++ir) {
      const double x = r[ir - 1];
      yy[ir][0] = std::pow(x, power) * (f0 + x * (f1 + x * f2));
      yy[ir][1] = std::pow(x, power) * (g0 + x * (g1 + x * g2));
    }
    return;
  }
  // Finite nucleus: potential regular at the origin, f, g ~ r^|k|.
  vzero = v1;
  if (kcur < 0) {
    ig = -kcur;
    f0 = 0.0;
    g0 = 1.0;
  } else {
    ig = kcur;
    f0 = 1.0;
    g0 = 0.0;
  }
  const double f1 = -(ecur - vzero) * abyt * g0 / double(ig - kcur + 1);
  const double g1 = (ecur - vzero + tbya * tbya) * abyt * f0 / double(ig + kcur + 1);
  const double f2 = -(ecur - vzero) * abyt * g1 / double(ig - kcur + 2);
  const double g2 = (ecur - vzero + tbya * tbya) * abyt * f1 / double(ig + kcur + 2);
  for (int ir = 1; ir <= 5; ++ir) {
    const double x = r[ir - 1];
    // integer power by repeated multiplication, as Fortran's r**ig
    const double xp = __builtin_powi(x, ig);
    yy[ir][0] = xp * (f0 + x * (f1 + x * f2));
    yy[ir][1] = xp * (g0 + x * (g1 + x * g2));
  }
}

RadialSolution DiracSolver::Solve(const atom::Orbital& orbital,
                                  std::span<const double> ruae0,
                                  double e0) const {
  const int mesh = static_cast<int>(grid_.size());
  const int ncur = orbital.n, lcur = orbital.l;
  const int kcur = orbital.Kappa();
  const double tbya = tbya_, abyt = abyt_;
  std::vector<double> r(mesh + 1), rab(mesh + 1), ruae(mesh + 1);
  for (int i = 1; i <= mesh; ++i) {
    r[i] = grid_.r(i - 1);
    rab[i] = grid_.rab()[i - 1];
    ruae[i] = ruae0[i - 1];
  }

  RadialSolution out;
  out.energy = e0;
  out.spinor = Spinor{Field(mesh, 0.0), Field(mesh, 0.0)};

  Coefficients zz;
  zz.z.resize(mesh + 1);
  for (int ir = 1; ir <= mesh; ++ir) {
    zz.z[ir][0][0] = rab[ir] * kcur / r[ir];
    zz.z[ir][1][1] = -zz.z[ir][0][0];
  }
  std::vector<std::array<double, 2>> yy(mesh + 1);

  double emin = -1.0e10, emax = 1.0, ecur = e0;
  double factor = 1.0;
  int ninf = mesh;
  for (int iter = 1; iter <= kMaxIterations; ++iter) {
    for (auto& y : yy) y = {0.0, 0.0};
    for (int ir = 1; ir <= mesh; ++ir) {
      zz.z[ir][0][1] = -rab[ir] * (ecur - ruae[ir]) * abyt;
      zz.z[ir][1][0] = -zz.z[ir][0][1] + rab[ir] * tbya;
    }
    // Classical turning point and practical infinity.
    int nctp = mesh;
    for (; nctp >= 10; --nctp) {
      if (zz.z[nctp][0][1] < 0.0) break;
    }
    if (nctp < 10) throw Error("DiracSolver: no classical turning point found");
    if (nctp > mesh - 10) {
      out.energy = e0 - 0.2;
      return out;
    }
    const double tolinf = std::pow(std::log(kThreshold), 2);
    for (ninf = nctp + 10; ninf <= mesh; ++ninf) {
      const double alpha2 = (ruae[ninf] - ecur) * std::pow(r[ninf] - r[nctp], 2);
      if (alpha2 > tolinf) break;
    }
    if (ninf > mesh) ninf = mesh;

    int ig = 0;
    StartAtOrigin(kcur, ecur, ruae[1], yy, ig);

    // Outward to the turning point.
    PredictorCorrector(zz, yy, 6, nctp);
    const double gout = yy[nctp][1];
    const double gpout =
        (zz.z[nctp][1][0] * yy[nctp][0] + zz.z[nctp][1][1] * yy[nctp][1]) / rab[nctp];

    // Inward from practical infinity.
    for (int ir = ninf; ir >= ninf - 4; --ir) {
      const double alpha = std::sqrt(ruae[ir] - ecur);
      yy[ir][1] = std::exp(-alpha * (r[ir] - r[nctp]));
      yy[ir][0] = (double(kcur) / r[ir] - alpha) * yy[ir][1] * tbya /
                  (ecur - ruae[ir] + tbya * tbya);
    }
    PredictorCorrector(zz, yy, ninf - 5, nctp);
    const double gin = yy[nctp][1];
    double gpin =
        (zz.z[nctp][1][0] * yy[nctp][0] + zz.z[nctp][1][1] * yy[nctp][1]) / rab[nctp];
    const double match = gout / gin;
    for (int ir = nctp; ir <= ninf; ++ir) {
      yy[ir][0] *= match;
      yy[ir][1] *= match;
    }
    gpin *= match;

    int nodes = 0;
    for (int i = 2; i <= ninf; ++i) {
      if (yy[i - 1][1] * yy[i][1] < 0.0) ++nodes;
    }
    const int wanted = ncur - lcur - 1;
    bool retry = false;
    if (nodes < wanted) {
      emin = ecur;
      ecur = (ecur * 0.9 > emax) ? 0.5 * ecur + 0.5 * emax : 0.9 * ecur;
      retry = true;
    } else if (nodes > wanted) {
      emax = ecur;
      ecur = (ecur * 1.1 < emin) ? 0.5 * ecur + 0.5 * emin : 1.1 * ecur;
      retry = true;
    }
    if (!retry) {
      Field dens(ninf);
      for (int ir = 1; ir <= ninf; ++ir) {
        dens[ir - 1] = yy[ir][0] * yy[ir][0] + yy[ir][1] * yy[ir][1];
      }
      factor = grid_.Integrate(dens, 2 * ig, ninf);
      const double decur = gout * (gpout - gpin) / factor;
      double decurp;
      if (decur > 0.0) {
        emin = ecur;
        decurp = std::min({decur, -0.2 * ecur, 0.7 * (emax - ecur)});
      } else {
        emax = ecur;
        decurp = -std::min({-decur, -0.2 * ecur, 0.7 * (ecur - emin)});
      }
      if (std::abs(decur) < kThreshold) {
        const double norm = std::sqrt(factor);
        for (int ir = 1; ir <= mesh; ++ir) {
          out.spinor.large[ir - 1] = yy[ir][1] / norm;
          out.spinor.small[ir - 1] = yy[ir][0] / norm;
        }
        out.energy = ecur;
        out.converged = true;
        return out;
      }
      ecur += decurp;
    }
  }
  out.energy = ecur;
  return out;
}

}  // namespace

std::unique_ptr<RadialSolver> MakeDiracSolver(const grid::RadialGrid& grid,
                                              const atom::Nucleus& nucleus,
                                              double speed_of_light) {
  return std::make_unique<DiracSolver>(grid, nucleus, speed_of_light);
}

}  // namespace atomx::solvers
