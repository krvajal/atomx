#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "atomx/atom/atom_state.hpp"
#include "atomx/core/types.hpp"

namespace atomx::xc {

// What an exchange-correlation model contributes for the current orbitals.
struct XcContribution {
  SpinField potential;  // Ry, per spin
  double energy = 0.0;  // Ry
  // Constant to subtract from printed eigenvalues of each spin (used by the
  // OEP, whose potential is only defined up to a constant).
  std::array<double, 2> eigenvalue_shift{0.0, 0.0};
};

// An exchange and/or correlation functional. Implementations may keep state
// between SCF iterations (e.g. the OEP warm start), hence non-const Evaluate.
class XcModel {
 public:
  virtual ~XcModel() = default;
  virtual std::string Name() const = 0;
  // `rho` is the radial density 4 pi r^2 n_s(r) of each spin.
  virtual XcContribution Evaluate(const atom::AtomState& state,
                                  const SpinField& rho) = 0;
};

// Sum of several models, e.g. KLI exchange + PZ correlation.
class CompositeXcModel final : public XcModel {
 public:
  CompositeXcModel& Add(std::unique_ptr<XcModel> model) {
    parts_.push_back(std::move(model));
    return *this;
  }
  std::string Name() const override;
  XcContribution Evaluate(const atom::AtomState& state,
                          const SpinField& rho) override;

 private:
  std::vector<std::unique_ptr<XcModel>> parts_;
};

}  // namespace atomx::xc
