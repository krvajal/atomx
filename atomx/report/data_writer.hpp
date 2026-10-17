#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atomx/calc/calculation.hpp"
#include "atomx/scf/scf_driver.hpp"

namespace atomx::report {

// Writes the radial data of a finished calculation (orbitals, densities,
// potentials, energies, SCF history) for plotting and post-processing.
class ResultWriter {
 public:
  virtual ~ResultWriter() = default;
  // Returns the paths of the files written.
  virtual std::vector<std::string> Write(const calc::Calculation& calculation,
                                         const scf::ScfResult& result) const = 0;
};

// Whitespace-separated columns with '#' headers, for gnuplot / xmgrace:
//   <prefix>.wfc    r, then r*R(r) of each orbital (large component)
//   <prefix>.wfcs   r, then the small components (Dirac only)
//   <prefix>.rho    r, 4 pi r^2 n_s(r) per spin
//   <prefix>.pot    r, V_nuc, V_H, V_xc,s, V_KS,s (Ry)
//   <prefix>.scf    iteration, residual
std::unique_ptr<ResultWriter> MakeColumnWriter(std::string prefix);

// Everything in one <prefix>.json, for Python / Jupyter.
std::unique_ptr<ResultWriter> MakeJsonWriter(std::string prefix);

// Writers named in a '+'- or ','-separated list: "dat", "json"; "none" or
// empty gives no writer.
std::vector<std::unique_ptr<ResultWriter>> MakeWriters(const std::string& list,
                                                       const std::string& prefix);

}  // namespace atomx::report
