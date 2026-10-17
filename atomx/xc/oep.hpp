#pragma once

#include "atomx/linalg/linear_solver.hpp"
#include "atomx/solvers/green_function.hpp"
#include "atomx/xc/exchange_operator.hpp"
#include "atomx/xc/xc_model.hpp"

namespace atomx::xc {

// Exact exchange-only optimized effective potential for non-relativistic
// atoms (ld1.x's dfx_new): a Krylov solve, in the Hartree metric, of the OEP
// equation \int chi_0(r,r') [Vx(r') - V_S(r')] = drho_0(r), preconditioned
// with the Thomas-Fermi response. Each call is warm-started from the previous
// OEP correction on top of the Slater potential of the current orbitals.
class OepExchange final : public XcModel {
 public:
  OepExchange(const ExchangeOperator& exchange,
              const solvers::OrbitalResponse& response,
              const linalg::LinearSolver& solver)
      : exchange_(exchange), response_(response), solver_(solver) {}

  std::string Name() const override { return "OEP"; }
  XcContribution Evaluate(const atom::AtomState& state,
                          const SpinField& rho) override;

 private:
  // Density response to a perturbation of each spin's potential.
  SpinField DensityResponse(const atom::AtomState& state,
                            const SpinField& dv) const;
  // Density response to the orbital-specific exchange terms K psi_k.
  SpinField ExchangeDensityResponse(const atom::AtomState& state,
                                    const std::vector<Spinor>& k_psi) const;

  const ExchangeOperator& exchange_;
  const solvers::OrbitalResponse& response_;
  const linalg::LinearSolver& solver_;

  // Warm start: OEP minus Slater from the previous call.
  SpinField previous_correction_;
  bool has_previous_ = false;
};

}  // namespace atomx::xc
