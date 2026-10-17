#include "atomx/scf/scf_driver.hpp"

#include <algorithm>
#include <cmath>

#include "atomx/core/constants.hpp"

namespace atomx::scf {

using constants::kE2;

SpinField ScfDriver::StartingPotential(atom::AtomState& state) const {
  const auto& orbitals = state.orbitals();
  const double z = nucleus_.charge();
  double enne = 0.0;
  for (std::size_t k = 0; k < state.size(); ++k) {
    enne += std::max(0.0, orbitals[k].occupation);
    double zen = 0.0;
    for (const auto& p : orbitals) {
      const double occ = std::max(0.0, p.occupation);
      if (p.n < orbitals[k].n) zen += occ;
      if (p.n == orbitals[k].n && p.l <= orbitals[k].l) zen += occ;
    }
    zen = std::max(z - zen + 1.0, 1.0);
    state.energy(k) = -std::pow(zen / orbitals[k].n, 2);
  }
  SpinField v{grid_.Zeros(), grid_.Zeros()};
  for (std::size_t i = 0; i < grid_.size(); ++i) {
    const double x = grid_.r(i) * std::cbrt(enne) / 0.885;
    double t = z / (1.0 + std::sqrt(x) * (0.02747 - x * (0.1486 - 0.007298 * x)) +
                    x * (1.243 + x * (0.2302 + 0.006944 * x)));
    t = std::max(1.0, t);
    v[0][i] = -kE2 * t / grid_.r(i);
  }
  v[1] = v[0];
  return v;
}

ScfResult ScfDriver::Run(const atom::Configuration& configuration,
                         Relativity relativity) {
  atom::AtomState state(grid_, configuration, relativity);
  const int nspin = state.nspin();
  const std::size_t mesh = grid_.size();
  const Field v_nucleus = nucleus_.Potential(grid_);
  SpinField vpot = StartingPotential(state);

  // Pack the potentials of the active spins into one vector for the mixer.
  auto pack = [&](const SpinField& v) {
    std::vector<double> p;
    p.reserve(mesh * nspin);
    for (int s = 0; s < nspin; ++s) p.insert(p.end(), v[s].begin(), v[s].end());
    return p;
  };
  auto unpack = [&](const std::vector<double>& p, SpinField& v) {
    for (int s = 0; s < nspin; ++s) {
      std::copy(p.begin() + s * mesh, p.begin() + (s + 1) * mesh, v[s].begin());
    }
  };

  int iterations = 0, failed = 0;
  double residual = 0.0;
  bool converged = false;
  Field v_hartree;
  SpinField rho;
  xc::XcContribution xc;
  for (int iter = 1; iter <= options_.max_iterations; ++iter) {
    // Orbitals in the current input potential.
    state.set_potential(vpot);
    int failures = 0;
    for (std::size_t k = 0; k < state.size(); ++k) {
      const auto& o = state.orbitals()[k];
      if (!o.IsActive()) {
        state.energy(k) = 0.0;
        state.spinor(k) = Spinor{grid_.Zeros(), grid_.Zeros()};
        continue;
      }
      auto sol = radial_solver_.Solve(o, vpot[o.spin], state.energy(k));
      if (!sol.converged) ++failures;
      state.energy(k) = sol.energy;
      state.spinor(k) = std::move(sol.spinor);
    }

    // New output potential from the new density.
    rho = state.Density();
    Field rho_total = rho[0];
    if (nspin == 2) {
      for (std::size_t i = 0; i < mesh; ++i) rho_total[i] += rho[1][i];
    }
    v_hartree = poisson_.Solve(0, 2, rho_total);
    for (double& v : v_hartree) v *= kE2;
    xc = xc_model_.Evaluate(state, rho);
    SpinField vnew{grid_.Zeros(), grid_.Zeros()};
    for (int s = 0; s < nspin; ++s) {
      for (std::size_t i = 0; i < mesh; ++i) {
        vnew[s][i] = v_nucleus[i] + v_hartree[i] + xc.potential[s][i];
      }
    }

    auto in = pack(vpot);
    const MixResult mix = mixer_->Mix(pack(vnew), in, iter);
    unpack(in, vpot);

    iterations = iter;
    residual = mix.residual;
    failed = failures;
    if (observer_) observer_(iter, mix.residual);
    if (mix.converged) {
      converged = true;
      break;
    }
  }

  const EnergyTerms energies =
      ComputeEnergies(state, rho, v_nucleus, v_hartree, xc);
  return ScfResult{std::move(state), std::move(rho), energies, std::move(xc),
                   iterations, residual, converged, failed};
}

}  // namespace atomx::scf
