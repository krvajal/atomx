#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atomx/linalg/dense_matrix.hpp"

namespace atomx::linalg {

enum class MatrixStructure {
  kGeneral,    // LU with partial pivoting
  kSymmetric,  // symmetric indefinite (Bunch-Kaufman or Aasen)
};

// Dense linear solves A x = b. The physics code only talks to this interface,
// so the backend (LAPACK on the host, SLATE on CPUs/GPUs/many nodes) can be
// swapped without touching it.
class LinearSolver {
 public:
  virtual ~LinearSolver() = default;

  virtual std::string name() const = 0;

  // Returns x. Throws atomx::Error if A is singular.
  virtual std::vector<double> Solve(const DenseMatrix& a,
                                    const std::vector<double>& b,
                                    MatrixStructure structure) const = 0;
};

// Host LAPACK (dgesv / dsysv).
std::unique_ptr<LinearSolver> MakeLapackSolver();

// Creates the solver selected by name: "lapack" or "slate" (when built with
// SLATE support).
std::unique_ptr<LinearSolver> MakeLinearSolver(const std::string& name);

}  // namespace atomx::linalg
