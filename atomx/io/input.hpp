#pragma once

#include <string>

#include "atomx/atom/configuration.hpp"
#include "atomx/core/types.hpp"
#include "atomx/io/namelist.hpp"

namespace atomx::io {

// Everything that defines an all-electron calculation. Field names and
// defaults follow ld1.x's &input namelist where an equivalent exists.
struct CalculationInput {
  std::string title;
  int z = 0;
  std::string configuration = "default";

  Relativity relativity = Relativity::kNone;
  bool spin_polarized = false;
  atom::SpinOrbitDistribution distribution = atom::SpinOrbitDistribution::kEnergy;

  // Exchange-correlation: "PZ" (= "SLA+PZ", LDA), "KLI", "OEP", or a '+'-joined
  // combination such as "KLI+PZ".
  std::string dft = "PZ";

  double xmin = -7.0;
  double dx = 0.008;
  double rmax = 100.0;

  double beta = 0.2;    // mixing
  double tr2 = 1e-14;   // SCF threshold on the mean squared potential change
  int max_iterations = 200;

  std::string nucleus = "point";  // or "sphere"
  double nuclear_radius = 0.0;    // bohr; 0 = default for the element
  double speed_of_light = 0.0;    // a.u.; 0 = CODATA value

  std::string linear_solver = "lapack";  // or "slate"

  // Reads &input (ld1.x names: zed/atom, config, rel, lsd, dft, xmin, dx, rmax,
  // beta, tr2, nucleus, rnucleus, cau_fact, rel_dist, title) plus the atomx
  // extension `linear_solver`. Throws on unsupported settings.
  static CalculationInput FromNamelist(const Namelist& nl);
};

}  // namespace atomx::io
