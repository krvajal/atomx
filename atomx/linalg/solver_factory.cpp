#include "atomx/core/types.hpp"
#include "atomx/linalg/linear_solver.hpp"

#ifdef ATOMX_WITH_SLATE
#include "atomx/linalg/slate_solver.hpp"
#endif

namespace atomx::linalg {

std::unique_ptr<LinearSolver> MakeLinearSolver(const std::string& name) {
  if (name == "lapack") return MakeLapackSolver();
#ifdef ATOMX_WITH_SLATE
  if (name == "slate") return MakeSlateSolver(SlateTarget::kHost);
  if (name == "slate-gpu") return MakeSlateSolver(SlateTarget::kDevices);
#else
  if (name == "slate" || name == "slate-gpu") {
    throw Error("atomx was built without SLATE; rebuild with --define slate=1");
  }
#endif
  throw Error("unknown linear solver '" + name + "'");
}

}  // namespace atomx::linalg
