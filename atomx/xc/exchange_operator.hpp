#pragma once

#include <memory>
#include <vector>

#include "atomx/atom/atom_state.hpp"
#include "atomx/core/types.hpp"
#include "atomx/grid/poisson.hpp"

namespace atomx::xc {

// The Fock exchange operator of the occupied orbitals, spherically averaged,
// applied to one orbital: returns [K psi_k](r) (Ry, same spinor layout as
// psi_k), defined so that E_x = 1/2 \sum_k occ_k <psi_k | K psi_k>.
class ExchangeOperator {
 public:
  virtual ~ExchangeOperator() = default;
  virtual Spinor Apply(const atom::AtomState& state, std::size_t k) const = 0;

  // K psi_k for every orbital, in orbital order.
  std::vector<Spinor> ApplyAll(const atom::AtomState& state) const;
  // E_x from precomputed K psi_k.
  static double Energy(const atom::AtomState& state,
                       const std::vector<Spinor>& k_psi);
  // <psi_k | K psi_k>
  static double Expectation(const atom::AtomState& state, std::size_t k,
                            const Spinor& k_psi);
};

// Non-relativistic exchange with ld1.x's spherical average for open shells
// (dvex): couples orbitals of the same spin through Legendre integrals.
std::unique_ptr<ExchangeOperator> MakeNonRelativisticExchange(
    const grid::RadialPoisson& poisson);

// Coulomb (longitudinal) exchange between Dirac spinors with j-coupled
// angular factors, Engel et al. PRA 58, 964 (1998), Eqs. (A15), (A17).
std::unique_ptr<ExchangeOperator> MakeDiracExchange(
    const grid::RadialPoisson& poisson);

}  // namespace atomx::xc
