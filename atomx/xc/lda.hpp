#pragma once

#include "atomx/xc/xc_model.hpp"

namespace atomx::xc {

// Local density functionals, as in Quantum ESPRESSO's XClib. Spin-polarised
// when the state has two spin channels.

// Slater exchange (alpha = 2/3).
class SlaterExchange final : public XcModel {
 public:
  std::string Name() const override { return "SLA"; }
  XcContribution Evaluate(const atom::AtomState& state,
                          const SpinField& rho) override;
};

// Perdew-Zunger (1981) correlation.
class PerdewZungerCorrelation final : public XcModel {
 public:
  std::string Name() const override { return "PZ"; }
  XcContribution Evaluate(const atom::AtomState& state,
                          const SpinField& rho) override;
};

}  // namespace atomx::xc
