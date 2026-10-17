#include "atomx/scf/energy.hpp"

namespace atomx::scf {

EnergyTerms ComputeEnergies(const atom::AtomState& state, const SpinField& rho,
                            const Field& v_nucleus, const Field& v_hartree,
                            const xc::XcContribution& xc) {
  const auto& grid = state.grid();
  const std::size_t mesh = grid.size();
  Field f_nuc(mesh), f_h(mesh), f_kin(mesh);
  for (std::size_t i = 0; i < mesh; ++i) {
    double rho_total = rho[0][i];
    if (state.nspin() == 2) rho_total += rho[1][i];
    f_nuc[i] = v_nucleus[i] * rho_total;
    f_h[i] = v_hartree[i] * rho_total;
    f_kin[i] = -f_nuc[i] - f_h[i];
    for (int s = 0; s < state.nspin(); ++s) {
      f_kin[i] -= xc.potential[s][i] * rho[s][i];
    }
  }
  EnergyTerms e;
  e.nuclear = grid.Integrate(f_nuc, 1);
  e.hartree = 0.5 * grid.Integrate(f_h, 2);
  e.xc = xc.energy;
  e.kinetic = grid.Integrate(f_kin, 1);
  for (std::size_t k = 0; k < state.size(); ++k) {
    const double occ = state.orbitals()[k].occupation;
    if (occ > 0.0) e.kinetic += occ * state.energy(k);
  }
  e.total = e.kinetic + e.nuclear + e.hartree + e.xc;
  return e;
}

}  // namespace atomx::scf
