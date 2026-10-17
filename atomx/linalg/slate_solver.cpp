#include "atomx/linalg/slate_solver.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <slate/slate.hh>
#include <string>
#include <vector>

#include "atomx/core/types.hpp"

namespace atomx::linalg {

namespace {

// Initialises MPI on first use and finalises it at program exit, unless the
// application manages MPI itself.
class MpiSession {
 public:
  static MpiSession& Get() {
    static MpiSession session;
    return session;
  }
  MPI_Comm comm() const { return MPI_COMM_WORLD; }
  int size() const { return size_; }

 private:
  MpiSession() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
      int provided = 0;
      MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
      if (provided < MPI_THREAD_MULTIPLE) {
        throw Error("SLATE needs MPI_THREAD_MULTIPLE");
      }
      owns_ = true;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
  }
  ~MpiSession() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (owns_ && !finalized) MPI_Finalize();
  }

  bool owns_ = false;
  int size_ = 1;
};

class SlateSolver final : public LinearSolver {
 public:
  explicit SlateSolver(SlateTarget target) : target_(target) {}

  std::string name() const override {
    return target_ == SlateTarget::kDevices ? "slate-gpu" : "slate";
  }

  std::vector<double> Solve(const DenseMatrix& a, const std::vector<double>& b,
                            MatrixStructure /*structure*/) const override {
    const int64_t n = static_cast<int64_t>(a.rows());
    if (a.cols() != a.rows() || b.size() != a.rows()) {
      throw Error("SlateSolver: dimension mismatch");
    }
    if (n == 0) return {};
    auto& mpi = MpiSession::Get();
    // The systems here are small (n < 100); one tile row/column per rank at
    // most. Larger problems would pick nb ~ 256-512.
    const int64_t nb = std::max<int64_t>(1, (n + mpi.size() - 1) / mpi.size());
    const int64_t p = mpi.size(), q = 1;
    const auto target = target_ == SlateTarget::kDevices
                            ? slate::Target::Devices
                            : slate::Target::HostTask;

    slate::Matrix<double> A(n, n, nb, p, q, mpi.comm());
    slate::Matrix<double> B(n, 1, nb, p, q, mpi.comm());
    A.insertLocalTiles();
    B.insertLocalTiles();
    for (int64_t j = 0; j < A.nt(); ++j) {
      for (int64_t i = 0; i < A.mt(); ++i) {
        if (!A.tileIsLocal(i, j)) continue;
        auto t = A(i, j);
        for (int64_t jj = 0; jj < t.nb(); ++jj) {
          for (int64_t ii = 0; ii < t.mb(); ++ii) {
            t.at(ii, jj) = a(i * nb + ii, j * nb + jj);
          }
        }
      }
    }
    for (int64_t i = 0; i < B.mt(); ++i) {
      if (!B.tileIsLocal(i, 0)) continue;
      auto t = B(i, 0);
      for (int64_t ii = 0; ii < t.mb(); ++ii) t.at(ii, 0) = b[i * nb + ii];
    }

    // LU with partial pivoting for both general and symmetric systems: the
    // OEP's symmetric matrix can be near-singular, where pivoted LU is the
    // safest choice.
    try {
      slate::lu_solve(A, B, {{slate::Option::Target, target}});
    } catch (const std::exception& e) {
      throw Error(std::string("SlateSolver: ") + e.what());
    }

    // Gather the distributed solution on every rank.
    std::vector<double> x(n, 0.0);
    for (int64_t i = 0; i < B.mt(); ++i) {
      if (!B.tileIsLocal(i, 0)) continue;
      B.tileGetForReading(i, 0, slate::LayoutConvert::ColMajor);
      auto t = B(i, 0);
      for (int64_t ii = 0; ii < t.mb(); ++ii) x[i * nb + ii] = t.at(ii, 0);
    }
    MPI_Allreduce(MPI_IN_PLACE, x.data(), static_cast<int>(n), MPI_DOUBLE,
                  MPI_SUM, mpi.comm());
    return x;
  }

 private:
  SlateTarget target_;
};

}  // namespace

std::unique_ptr<LinearSolver> MakeSlateSolver(SlateTarget target) {
  return std::make_unique<SlateSolver>(target);
}

}  // namespace atomx::linalg
