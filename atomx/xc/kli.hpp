#pragma once

#include "atomx/linalg/linear_solver.hpp"
#include "atomx/xc/exchange_operator.hpp"
#include "atomx/xc/xc_model.hpp"

namespace atomx::xc {

// Krieger-Li-Iafrate approximation to the exchange-only OEP,
//   Vx(r) = V_S(r) + \sum_{i != HOMO} occ_i n_i(r) / rho(r) x_i,
// with the constants x_i from a small linear system. Works with any exchange
// operator, so the same class gives non-relativistic KLI and relativistic
// KLI (Kreibich, Gross, Engel, PRA 57, 138 (1998), Eq. 54).
class KliExchange final : public XcModel {
 public:
  KliExchange(const ExchangeOperator& exchange,
              const linalg::LinearSolver& solver)
      : exchange_(exchange), solver_(solver) {}

  std::string Name() const override { return "KLI"; }
  XcContribution Evaluate(const atom::AtomState& state,
                          const SpinField& rho) override;

 private:
  const ExchangeOperator& exchange_;
  const linalg::LinearSolver& solver_;
};

}  // namespace atomx::xc
