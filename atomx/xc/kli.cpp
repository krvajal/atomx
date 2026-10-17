#include "atomx/xc/kli.hpp"

#include <limits>
#include <vector>

namespace atomx::xc {

namespace {

// Last grid point of the region, contiguous from the origin, where rho > 0.
// The radial solvers set far tails exactly to zero, so 1/rho is only formed
// below this point.
std::size_t LastPositive(const Field& rho) {
  std::size_t n = 0;
  while (n < rho.size() && rho[n] > std::numeric_limits<double>::min()) ++n;
  return n;
}

}  // namespace

XcContribution KliExchange::Evaluate(const atom::AtomState& state,
                                     const SpinField& rho) {
  const auto& grid = state.grid();
  const auto& orbitals = state.orbitals();
  const std::size_t mesh = grid.size();
  const auto k_psi = exchange_.ApplyAll(state);

  XcContribution out;
  out.potential = {grid.Zeros(), grid.Zeros()};
  out.energy = ExchangeOperator::Energy(state, k_psi);

  for (int s = 0; s < state.nspin(); ++s) {
    // Occupied orbitals of this spin; the HOMO has the highest eigenvalue.
    std::vector<std::size_t> occ;
    std::size_t homo = 0;
    for (std::size_t k = 0; k < state.size(); ++k) {
      if (orbitals[k].spin != s || !orbitals[k].IsOccupied()) continue;
      if (occ.empty() || state.energy(k) > state.energy(occ[homo])) {
        homo = occ.size();
      }
      occ.push_back(k);
    }
    if (occ.empty()) continue;
    const std::size_t last = LastPositive(rho[s]);

    std::vector<Field> dens(occ.size());
    Field slater(mesh, 0.0);
    for (std::size_t a = 0; a < occ.size(); ++a) {
      const std::size_t k = occ[a];
      dens[a] = state.OrbitalDensity(k);
      const auto& psi = state.spinor(k);
      const double w = orbitals[k].occupation;
      for (std::size_t i = 0; i < last; ++i) {
        slater[i] += w * (psi.large[i] * k_psi[k].large[i] +
                          psi.small[i] * k_psi[k].small[i]);
      }
    }
    for (std::size_t i = 0; i < last; ++i) slater[i] /= rho[s][i];

    // <i|u_i|i> and <i|V_S|i>
    std::vector<double> u_avg(occ.size()), vs_avg(occ.size());
    for (std::size_t a = 0; a < occ.size(); ++a) {
      const std::size_t k = occ[a];
      const auto& psi = state.spinor(k);
      Field w(mesh);
      for (std::size_t i = 0; i < mesh; ++i) {
        w[i] = psi.large[i] * k_psi[k].large[i] +
               psi.small[i] * k_psi[k].small[i];
      }
      const int nst = state.DensityPower(k);
      u_avg[a] = grid.Integrate(w, nst, last);
      for (std::size_t i = 0; i < mesh; ++i) w[i] = dens[a][i] * slater[i];
      vs_avg[a] = grid.Integrate(w, nst, last);
    }

    // (I - M) x = <V_S> - <u> over the non-HOMO orbitals.
    std::vector<std::size_t> rows;
    for (std::size_t a = 0; a < occ.size(); ++a) {
      if (a != homo) rows.push_back(a);
    }
    std::vector<double> x_of(occ.size(), 0.0);
    if (!rows.empty()) {
      const std::size_t n = rows.size();
      linalg::DenseMatrix m(n, n);
      std::vector<double> y(n);
      Field w(mesh, 0.0);
      for (std::size_t p = 0; p < n; ++p) {
        const std::size_t aj = rows[p];
        const int nst = state.DensityPower(occ[aj]);
        for (std::size_t q = 0; q < n; ++q) {
          const std::size_t ai = rows[q];
          const double occ_i = orbitals[occ[ai]].occupation;
          std::fill(w.begin(), w.end(), 0.0);
          for (std::size_t i = 0; i < last; ++i) {
            w[i] = dens[aj][i] * occ_i * dens[ai][i] / rho[s][i];
          }
          m(p, q) = -grid.Integrate(w, nst, last);
        }
        m(p, p) += 1.0;
        y[p] = vs_avg[aj] - u_avg[aj];
      }
      const auto x =
          solver_.Solve(m, y, linalg::MatrixStructure::kGeneral);
      for (std::size_t p = 0; p < n; ++p) x_of[rows[p]] = x[p];
    }

    Field& v = out.potential[s];
    for (std::size_t i = 0; i < last; ++i) {
      double corr = 0.0;
      for (std::size_t a = 0; a < occ.size(); ++a) {
        corr += dens[a][i] * x_of[a] * orbitals[occ[a]].occupation;
      }
      v[i] = slater[i] + corr / rho[s][i];
    }
  }
  return out;
}

}  // namespace atomx::xc
