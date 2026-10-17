#pragma once

#include <ostream>

#include "atomx/calc/calculation.hpp"
#include "atomx/scf/scf_driver.hpp"

namespace atomx::report {

// Prints the results in the layout of ld1.x's output (energies in Ry, Ha, eV
// and the eigenvalue table), so the two programs can be compared with diff.
void WriteReport(std::ostream& out, const calc::Calculation& calculation,
                 const scf::ScfResult& result);

}  // namespace atomx::report
