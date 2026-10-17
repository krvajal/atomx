#pragma once

// Helpers shared by the tests.

#include <cmath>
#include <functional>
#include <string>

#include "atomx/atom/atom_state.hpp"
#include "atomx/atom/configuration.hpp"
#include "atomx/grid/radial_grid.hpp"
#include "atomx/io/input.hpp"
#include "atomx/io/namelist.hpp"

namespace atomx::testing {

// Fills orbital k of `state` with the radial function f(r) (as r*R(r), large
// component only), normalised on the grid.
inline void SetOrbital(atom::AtomState& state, std::size_t k,
                       const std::function<double(double)>& f) {
  const auto& g = state.grid();
  Field large(g.size()), dens(g.size());
  for (std::size_t i = 0; i < g.size(); ++i) {
    large[i] = f(g.r(i));
    dens[i] = large[i] * large[i];
  }
  const double norm = std::sqrt(g.Integrate(dens, state.DensityPower(k)));
  for (double& x : large) x /= norm;
  state.spinor(k) = Spinor{large, g.Zeros()};
}

// Hydrogenic radial shapes (unnormalised r*R) for nuclear charge z.
inline double Hydrogen1s(double z, double r) { return r * std::exp(-z * r); }
inline double Hydrogen2s(double z, double r) {
  return r * (1.0 - 0.5 * z * r) * std::exp(-0.5 * z * r);
}
inline double Hydrogen2p(double z, double r) {
  return r * r * std::exp(-0.5 * z * r);
}

// A CalculationInput from the body of an ld1.x-style &input namelist.
inline io::CalculationInput InputFromNamelist(const std::string& body) {
  return io::CalculationInput::FromNamelist(
      io::Namelist::FromString("&input\n" + body + "\n/\n", "input"));
}

}  // namespace atomx::testing
