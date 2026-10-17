#pragma once

#include <cstddef>
#include <vector>

namespace atomx::linalg {

// Small column-major dense matrix: the exchange of data between the physics
// code and a linear-algebra backend. Backends copy it into their own layout
// (e.g. a distributed slate::Matrix) as needed.
class DenseMatrix {
 public:
  DenseMatrix() = default;
  DenseMatrix(std::size_t rows, std::size_t cols)
      : rows_(rows), cols_(cols), data_(rows * cols, 0.0) {}

  std::size_t rows() const { return rows_; }
  std::size_t cols() const { return cols_; }

  double& operator()(std::size_t i, std::size_t j) {
    return data_[i + j * rows_];
  }
  double operator()(std::size_t i, std::size_t j) const {
    return data_[i + j * rows_];
  }

  double* data() { return data_.data(); }
  const double* data() const { return data_.data(); }

 private:
  std::size_t rows_ = 0;
  std::size_t cols_ = 0;
  std::vector<double> data_;
};

}  // namespace atomx::linalg
