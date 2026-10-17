#pragma once

#include <memory>

#include "atomx/atom/configuration.hpp"
#include "atomx/atom/nucleus.hpp"
#include "atomx/grid/poisson.hpp"
#include "atomx/grid/radial_grid.hpp"
#include "atomx/io/input.hpp"
#include "atomx/linalg/linear_solver.hpp"
#include "atomx/scf/scf_driver.hpp"
#include "atomx/solvers/green_function.hpp"
#include "atomx/solvers/radial_solver.hpp"
#include "atomx/xc/exchange_operator.hpp"
#include "atomx/xc/xc_model.hpp"

namespace atomx::calc {

// Composition root: builds every component from a CalculationInput and owns
// them for the lifetime of the run. Swapping a component (another solver,
// functional, nucleus or linear-algebra backend) only touches this class.
class Calculation {
 public:
  explicit Calculation(io::CalculationInput input);
  Calculation(const Calculation&) = delete;
  Calculation& operator=(const Calculation&) = delete;

  scf::ScfResult Run(scf::ScfDriver::Observer observer = {});

  const io::CalculationInput& input() const { return input_; }
  const grid::RadialGrid& grid() const { return grid_; }
  const atom::Nucleus& nucleus() const { return *nucleus_; }
  const atom::Configuration& configuration() const { return configuration_; }
  const xc::XcModel& xc_model() const { return *xc_model_; }
  double speed_of_light() const { return speed_of_light_; }

 private:
  std::unique_ptr<xc::XcModel> BuildXcModel(const std::string& dft);

  io::CalculationInput input_;
  grid::RadialGrid grid_;
  grid::RadialPoisson poisson_;
  std::unique_ptr<atom::Nucleus> nucleus_;
  atom::Configuration configuration_;
  double speed_of_light_;
  std::unique_ptr<linalg::LinearSolver> linear_solver_;
  std::unique_ptr<xc::ExchangeOperator> exchange_;
  std::unique_ptr<solvers::RadialSolver> radial_solver_;
  std::unique_ptr<solvers::OrbitalResponse> response_;
  std::unique_ptr<xc::XcModel> xc_model_;
};

}  // namespace atomx::calc
