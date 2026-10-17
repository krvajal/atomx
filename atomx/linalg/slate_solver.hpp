#pragma once

#include <memory>

#include "atomx/linalg/linear_solver.hpp"

namespace atomx::linalg {

enum class SlateTarget {
  kHost,     // multithreaded CPU tasks
  kDevices,  // GPUs (CUDA), tiles moved to the devices by SLATE
};

// Dense solves with SLATE (https://github.com/icl-utk-edu/slate) on a
// p x 1 process grid over MPI_COMM_WORLD. Initialises MPI if the caller has
// not. Every rank must call Solve with the same data, in the same order
// (the atomx SCF is deterministic, so each rank runs it redundantly).
std::unique_ptr<LinearSolver> MakeSlateSolver(SlateTarget target);

}  // namespace atomx::linalg
