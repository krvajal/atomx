#pragma once

#include "atomx/atom/atom_state.hpp"
#include "atomx/core/types.hpp"
#include "atomx/xc/xc_model.hpp"

namespace atomx::scf {

// Total-energy decomposition, in Ry (ld1.x's elsd).
struct EnergyTerms {
  double kinetic = 0.0;    // Ekin = sum occ e_k - \int V_KS rho
  double nuclear = 0.0;    // Encl = \int V_nuc rho
  double hartree = 0.0;    // Eh
  double xc = 0.0;         // Exc
  double total = 0.0;

  // Levy-Perdew virial discrepancy -(Ekin + Etot), zero for exact OEP.
  double VirialDiscrepancy() const { return -(kinetic + total); }
};

EnergyTerms ComputeEnergies(const atom::AtomState& state, const SpinField& rho,
                            const Field& v_nucleus, const Field& v_hartree,
                            const xc::XcContribution& xc);

}  // namespace atomx::scf
