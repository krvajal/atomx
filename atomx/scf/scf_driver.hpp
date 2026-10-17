#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "atomx/atom/atom_state.hpp"
#include "atomx/atom/configuration.hpp"
#include "atomx/atom/nucleus.hpp"
#include "atomx/grid/poisson.hpp"
#include "atomx/scf/energy.hpp"
#include "atomx/scf/mixer.hpp"
#include "atomx/solvers/radial_solver.hpp"
#include "atomx/xc/xc_model.hpp"

namespace atomx::scf {

struct ScfOptions {
  int max_iterations = 200;
};

struct ScfResult {
  atom::AtomState state;
  SpinField density;
  EnergyTerms energies;
  xc::XcContribution xc;
  int iterations = 0;
  double residual = 0.0;
  bool converged = false;
  int failed_orbital_solves = 0;  // in the last iteration
  // Potentials of the last iteration (Ry) and the residual of every step.
  Field v_nucleus;
  Field v_hartree;
  std::vector<double> residual_history;
};

// Kohn-Sham self-consistency for a spherical atom. All physics is injected:
// the radial solver (Schroedinger or Dirac), the XC model, the mixer and the
// nucleus, so any combination can be run by the same loop.
class ScfDriver {
 public:
  ScfDriver(const grid::RadialGrid& grid, const atom::Nucleus& nucleus,
            const solvers::RadialSolver& radial_solver, xc::XcModel& xc_model,
            std::unique_ptr<Mixer> mixer, ScfOptions options = {})
      : grid_(grid),
        poisson_(grid),
        nucleus_(nucleus),
        radial_solver_(radial_solver),
        xc_model_(xc_model),
        mixer_(std::move(mixer)),
        options_(options) {}

  // Called after every iteration (for progress output); may be empty.
  using Observer = std::function<void(int iteration, double residual)>;
  void SetObserver(Observer observer) { observer_ = std::move(observer); }

  ScfResult Run(const atom::Configuration& configuration,
                Relativity relativity);

 private:
  // Thomas-Fermi starting potential and hydrogenic eigenvalue guesses
  // (ld1.x's starting_potential).
  SpinField StartingPotential(atom::AtomState& state) const;

  const grid::RadialGrid& grid_;
  grid::RadialPoisson poisson_;
  const atom::Nucleus& nucleus_;
  const solvers::RadialSolver& radial_solver_;
  xc::XcModel& xc_model_;
  std::unique_ptr<Mixer> mixer_;
  ScfOptions options_;
  Observer observer_;
};

}  // namespace atomx::scf
