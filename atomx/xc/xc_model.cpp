#include "atomx/xc/xc_model.hpp"

namespace atomx::xc {

std::string CompositeXcModel::Name() const {
  std::string name;
  for (const auto& p : parts_) {
    if (!name.empty()) name += "+";
    name += p->Name();
  }
  return name.empty() ? "none" : name;
}

XcContribution CompositeXcModel::Evaluate(const atom::AtomState& state,
                                          const SpinField& rho) {
  const auto& grid = state.grid();
  XcContribution total;
  total.potential = {grid.Zeros(), grid.Zeros()};
  for (const auto& p : parts_) {
    const XcContribution c = p->Evaluate(state, rho);
    for (int s = 0; s < 2; ++s) {
      for (std::size_t i = 0; i < grid.size(); ++i) {
        total.potential[s][i] += c.potential[s][i];
      }
      total.eigenvalue_shift[s] += c.eigenvalue_shift[s];
    }
    total.energy += c.energy;
  }
  return total;
}

}  // namespace atomx::xc
