#pragma once

#include <cmath>
#include <vector>

#include "atomx/atom/configuration.hpp"
#include "atomx/atom/orbital.hpp"
#include "atomx/core/types.hpp"
#include "atomx/grid/radial_grid.hpp"

namespace atomx::atom {

// The Kohn-Sham orbitals of an atom at one point of the SCF cycle: the
// orbital list, their eigenvalues (Ry) and radial spinors r*R(r), normalised
// as \int (G^2 + F^2) dr = 1.
class AtomState {
 public:
  AtomState(const grid::RadialGrid& grid, const Configuration& configuration,
            Relativity relativity)
      : grid_(grid),
        orbitals_(configuration.orbitals()),
        nspin_(configuration.nspin()),
        relativity_(relativity),
        energies_(orbitals_.size(), 0.0),
        spinors_(orbitals_.size(),
                 Spinor{grid.Zeros(), grid.Zeros()}) {}

  const grid::RadialGrid& grid() const { return grid_; }
  const std::vector<Orbital>& orbitals() const { return orbitals_; }
  std::size_t size() const { return orbitals_.size(); }
  int nspin() const { return nspin_; }
  Relativity relativity() const { return relativity_; }

  double energy(std::size_t k) const { return energies_[k]; }
  double& energy(std::size_t k) { return energies_[k]; }
  const Spinor& spinor(std::size_t k) const { return spinors_[k]; }
  Spinor& spinor(std::size_t k) { return spinors_[k]; }

  // The effective (Kohn-Sham) potential, per spin, of which the current
  // orbitals are eigenstates. Needed by response-based models (OEP).
  const SpinField& potential() const { return potential_; }
  void set_potential(const SpinField& v) { potential_ = v; }

  // G^2 + F^2 of orbital k.
  Field OrbitalDensity(std::size_t k) const {
    Field d(grid_.size());
    const auto& s = spinors_[k];
    for (std::size_t i = 0; i < d.size(); ++i) {
      d[i] = s.large[i] * s.large[i] + s.small[i] * s.small[i];
    }
    return d;
  }

  // Radial density 4 pi r^2 n_s(r) of each spin channel.
  SpinField Density() const {
    SpinField rho{grid_.Zeros(), grid_.Zeros()};
    for (std::size_t k = 0; k < size(); ++k) {
      const auto& o = orbitals_[k];
      if (!o.IsActive()) continue;
      const auto& s = spinors_[k];
      auto& r = rho[o.spin];
      for (std::size_t i = 0; i < r.size(); ++i) {
        r[i] += o.occupation *
                (s.large[i] * s.large[i] + s.small[i] * s.small[i]);
      }
    }
    return rho;
  }

  // Leading power of an orbital density at small r, for grid integration.
  int DensityPower(std::size_t k) const {
    const auto& o = orbitals_[k];
    if (relativity_ == Relativity::kDirac) {
      return static_cast<int>(std::lround(2.0 * o.j + 1.0));
    }
    return 2 * o.l + 2;
  }

 private:
  const grid::RadialGrid& grid_;
  std::vector<Orbital> orbitals_;
  int nspin_;
  Relativity relativity_;
  std::vector<double> energies_;
  std::vector<Spinor> spinors_;
  SpinField potential_;
};

}  // namespace atomx::atom
