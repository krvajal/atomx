#include "atomx/xc/oep.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "atomx/core/constants.hpp"

namespace atomx::xc {

namespace {

constexpr int kMaxKrylov = 12;
// Converged when the residual norm (Hartree metric) is below this.
constexpr double kTolerance = 1.0e-11;

std::size_t LastPositive(const Field& rho) {
  std::size_t n = 0;
  while (n < rho.size() && rho[n] > std::numeric_limits<double>::min()) ++n;
  return n;
}

}  // namespace

SpinField OepExchange::DensityResponse(const atom::AtomState& state,
                                       const SpinField& dv) const {
  const auto& grid = state.grid();
  SpinField drho{grid.Zeros(), grid.Zeros()};
  for (std::size_t k = 0; k < state.size(); ++k) {
    const auto& o = state.orbitals()[k];
    if (!o.IsOccupied()) continue;
    const auto& psi = state.spinor(k).large;
    Field dvy(grid.size()), w(grid.size());
    for (std::size_t i = 0; i < dvy.size(); ++i) {
      dvy[i] = dv[o.spin][i] * psi[i];
      w[i] = dvy[i] * psi[i];
    }
    const double fac = grid.Integrate(w, 2 * o.l + 2);
    for (std::size_t i = 0; i < dvy.size(); ++i) dvy[i] -= fac * psi[i];
    const Field dchi = response_.Solve(o.l, state.energy(k), dvy, psi,
                                       state.potential()[o.spin]);
    for (std::size_t i = 0; i < dvy.size(); ++i) {
      drho[o.spin][i] += 2.0 * o.occupation * psi[i] * dchi[i];
    }
  }
  return drho;
}

SpinField OepExchange::ExchangeDensityResponse(
    const atom::AtomState& state, const std::vector<Spinor>& k_psi) const {
  const auto& grid = state.grid();
  SpinField drho{grid.Zeros(), grid.Zeros()};
  for (std::size_t k = 0; k < state.size(); ++k) {
    const auto& o = state.orbitals()[k];
    if (!o.IsOccupied()) continue;
    const auto& psi = state.spinor(k).large;
    Field dvy = k_psi[k].large, w(grid.size());
    for (std::size_t i = 0; i < w.size(); ++i) w[i] = dvy[i] * psi[i];
    const double fac = grid.Integrate(w, 2 * o.l + 2);
    for (std::size_t i = 0; i < dvy.size(); ++i) dvy[i] -= fac * psi[i];
    const Field dchi = response_.Solve(o.l, state.energy(k), dvy, psi,
                                       state.potential()[o.spin]);
    for (std::size_t i = 0; i < dvy.size(); ++i) {
      drho[o.spin][i] += 2.0 * o.occupation * psi[i] * dchi[i];
    }
  }
  return drho;
}

XcContribution OepExchange::Evaluate(const atom::AtomState& state,
                                     const SpinField& rho) {
  const auto& grid = state.grid();
  const auto& orbitals = state.orbitals();
  const std::size_t mesh = grid.size();
  const int nspin = state.nspin();
  const grid::RadialPoisson poisson(grid);
  const auto k_psi = exchange_.ApplyAll(state);

  // <x, V_H[y]> summed over spins (ld1.x weights spins by 2 when polarised).
  auto hartree = [&](const SpinField& f) {
    SpinField v{grid.Zeros(), grid.Zeros()};
    for (int s = 0; s < nspin; ++s) v[s] = poisson.Solve(0, 2, f[s]);
    return v;
  };
  auto metric = [&](const SpinField& x, const SpinField& vh) {
    Field aux(mesh);
    for (std::size_t i = 0; i < mesh; ++i) {
      aux[i] = x[0][i] * vh[0][i];
      if (nspin == 2) aux[i] = 2.0 * (aux[i] + x[1][i] * vh[1][i]);
    }
    return grid.Integrate(aux, 2);
  };

  std::array<std::size_t, 2> last{0, 0};
  for (int s = 0; s < nspin; ++s) last[s] = LastPositive(rho[s]);

  const SpinField drho0 = ExchangeDensityResponse(state, k_psi);

  // Slater potential of the current orbitals, continued as 1/r beyond the
  // region where rho > 0.
  SpinField vs0{grid.Zeros(), grid.Zeros()};
  for (int s = 0; s < nspin; ++s) {
    for (std::size_t k = 0; k < state.size(); ++k) {
      if (orbitals[k].spin != s) continue;
      const auto& psi = state.spinor(k).large;
      for (std::size_t i = 0; i < mesh; ++i) {
        vs0[s][i] += orbitals[k].occupation * psi[i] * k_psi[k].large[i];
      }
    }
    const std::size_t n = last[s];
    if (n == 0) {
      std::fill(vs0[s].begin(), vs0[s].end(), 0.0);
      continue;
    }
    for (std::size_t i = 0; i < n; ++i) vs0[s][i] /= rho[s][i];
    for (std::size_t i = n; i < mesh; ++i) {
      vs0[s][i] = vs0[s][n - 1] * grid.r(n - 1) / grid.r(i);
    }
  }

  SpinField vslater = vs0;
  if (has_previous_) {
    for (int s = 0; s < nspin; ++s) {
      for (std::size_t i = 0; i < mesh; ++i) {
        vslater[s][i] += previous_correction_[s][i];
      }
    }
  }

  SpinField drho1 = DensityResponse(state, vslater);
  for (int s = 0; s < nspin; ++s) {
    for (std::size_t i = 0; i < mesh; ++i) drho1[s][i] -= drho0[s][i];
  }
  const SpinField dvh1 = hartree(drho1);
  const double c1 = metric(drho1, dvh1);

  // Thomas-Fermi approximation to the inverse response.
  const double fac = -std::pow(0.75 * constants::kPi, 2.0 / 3.0);
  SpinField appchim1{grid.Zeros(), grid.Zeros()};
  for (int s = 0; s < nspin; ++s) {
    for (std::size_t i = 0; i < last[s]; ++i) {
      const double r = grid.r(i);
      appchim1[s][i] = fac / (r * std::cbrt(nspin * rho[s][i] * r));
    }
  }

  SpinField vx = vslater;
  auto finish = [&]() {
    previous_correction_ = {grid.Zeros(), grid.Zeros()};
    for (int s = 0; s < nspin; ++s) {
      for (std::size_t i = 0; i < mesh; ++i) {
        previous_correction_[s][i] = vx[s][i] - vs0[s][i];
      }
    }
    has_previous_ = true;

    XcContribution out;
    out.potential = vx;
    out.energy = ExchangeOperator::Energy(state, k_psi);
    // Shift printed eigenvalues so that the HOMO condition holds.
    for (int s = 0; s < nspin; ++s) {
      double ehomo = -std::numeric_limits<double>::max();
      for (std::size_t k = 0; k < state.size(); ++k) {
        const auto& o = orbitals[k];
        if (o.spin != s || !o.IsOccupied() || state.energy(k) <= ehomo) continue;
        ehomo = state.energy(k);
        const auto& psi = state.spinor(k).large;
        Field w(mesh);
        for (std::size_t i = 0; i < mesh; ++i) w[i] = vx[s][i] * psi[i] * psi[i];
        out.eigenvalue_shift[s] = grid.Integrate(w, 2 * o.l + 2) -
                                  ExchangeOperator::Expectation(state, k, k_psi[k]);
      }
    }
    return out;
  };

  if (c1 < kTolerance) return finish();

  std::vector<SpinField> vvx(kMaxKrylov), drhox(kMaxKrylov + 1),
      dvh(kMaxKrylov + 1);
  drhox[0] = drho1;
  linalg::DenseMatrix a(kMaxKrylov, kMaxKrylov);
  std::vector<double> b1(kMaxKrylov), b2(kMaxKrylov);

  for (int it = 0; it < kMaxKrylov; ++it) {
    // New search direction: preconditioned residual, normalised.
    vvx[it] = {grid.Zeros(), grid.Zeros()};
    for (int s = 0; s < nspin; ++s) {
      for (std::size_t i = 0; i < mesh; ++i) {
        vvx[it][s][i] = appchim1[s][i] * drhox[it][s][i];
      }
    }
    dvh[it] = hartree(vvx[it]);
    const double norm = metric(vvx[it], dvh[it]);
    if (norm > 0.0) {
      const double scale = 1.0 / std::sqrt(norm);
      for (int s = 0; s < nspin; ++s) {
        for (double& x : vvx[it][s]) x *= scale;
      }
    }
    drhox[it] = DensityResponse(state, vvx[it]);
    dvh[it] = hartree(drhox[it]);
    b1[it] = metric(drhox[it], dvh1);
    b2[it] = metric(drho1, dvh[it]);
    for (int j = 0; j <= it; ++j) {
      a(it, j) = metric(drhox[it], dvh[j]);
      a(j, it) = metric(drhox[j], dvh[it]);
    }
    const int n = it + 1;
    linalg::DenseMatrix am(n, n);
    std::vector<double> rhs(n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        a(i, j) = a(j, i) = 0.5 * (a(i, j) + a(j, i));
      }
    }
    for (int i = 0; i < n; ++i) {
      rhs[i] = -0.5 * (b1[i] + b2[i]);
      for (int j = 0; j < n; ++j) am(i, j) = a(i, j);
    }
    const auto x = solver_.Solve(am, rhs, linalg::MatrixStructure::kSymmetric);

    // Residual norm of the new potential.
    double capel = c1;
    for (int i = 0; i < n; ++i) {
      capel += x[i] * (b1[i] + b2[i]);
      for (int j = 0; j < n; ++j) capel += x[i] * a(i, j) * x[j];
    }
    vx = vslater;
    for (int j = 0; j < n; ++j) {
      for (int s = 0; s < nspin; ++s) {
        for (std::size_t i = 0; i < mesh; ++i) vx[s][i] += x[j] * vvx[j][s][i];
      }
    }
    if (n == kMaxKrylov || capel < kTolerance) {
      // At the numerical noise floor, keep only real improvements.
      if (capel > 0.5 * c1) vx = vslater;
      return finish();
    }
    drhox[it + 1] = DensityResponse(state, vx);
    for (int s = 0; s < nspin; ++s) {
      for (std::size_t i = 0; i < mesh; ++i) drhox[it + 1][s][i] -= drho0[s][i];
    }
  }
  return finish();  // not reached
}

}  // namespace atomx::xc
