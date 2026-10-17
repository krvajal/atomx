#include <string>
#include <vector>

#include "atomx/core/types.hpp"
#include "atomx/linalg/linear_solver.hpp"

extern "C" {
void dgesv_(const int* n, const int* nrhs, double* a, const int* lda,
            int* ipiv, double* b, const int* ldb, int* info);
void dsysv_(const char* uplo, const int* n, const int* nrhs, double* a,
            const int* lda, int* ipiv, double* b, const int* ldb,
            double* work, const int* lwork, int* info);
}

namespace atomx::linalg {

namespace {

class LapackSolver final : public LinearSolver {
 public:
  std::string name() const override { return "lapack"; }

  std::vector<double> Solve(const DenseMatrix& a, const std::vector<double>& b,
                            MatrixStructure structure) const override {
    const int n = static_cast<int>(a.rows());
    if (a.cols() != a.rows() || b.size() != a.rows()) {
      throw Error("LapackSolver: dimension mismatch");
    }
    if (n == 0) return {};
    DenseMatrix lu = a;
    std::vector<double> x = b;
    std::vector<int> ipiv(n);
    const int nrhs = 1;
    int info = 0;
    if (structure == MatrixStructure::kGeneral) {
      dgesv_(&n, &nrhs, lu.data(), &n, ipiv.data(), x.data(), &n, &info);
    } else {
      const char uplo = 'U';
      const int lwork = 64 * n;
      std::vector<double> work(lwork);
      dsysv_(&uplo, &n, &nrhs, lu.data(), &n, ipiv.data(), x.data(), &n,
             work.data(), &lwork, &info);
    }
    if (info != 0) {
      throw Error("LapackSolver: singular matrix (info=" +
                  std::to_string(info) + ")");
    }
    return x;
  }
};

}  // namespace

std::unique_ptr<LinearSolver> MakeLapackSolver() {
  return std::make_unique<LapackSolver>();
}

}  // namespace atomx::linalg
