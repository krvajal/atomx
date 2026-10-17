#include "atomx/scf/mixer.hpp"

#include <cmath>
#include <numeric>

namespace atomx::scf {

namespace {

double Dot(const std::vector<double>& a, const std::vector<double>& b) {
  return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

class AndersonMixer final : public Mixer {
 public:
  AndersonMixer(double beta, double threshold)
      : beta_(beta), threshold_(threshold) {}

  MixResult Mix(const std::vector<double>& output, std::vector<double>& b,
                int iteration) override {
    const std::size_t n = b.size();
    std::vector<double> a(n);
    for (std::size_t i = 0; i < n; ++i) a[i] = output[i] - b[i];
    MixResult res;
    res.residual = Dot(a, a) / static_cast<double>(n);
    if (res.residual < threshold_) {
      res.converged = true;
      return res;
    }
    if (iteration == 1) {
      a1_ = a;
      b1_ = b;
      for (std::size_t i = 0; i < n; ++i) b[i] += beta_ * a[i];
      return res;
    }
    std::vector<double> c(n), d;
    for (std::size_t i = 0; i < n; ++i) c[i] = a1_[i] - a[i];
    const double d11 = Dot(c, c);
    const double rd1m = Dot(a, c);
    std::vector<double> next(n);
    if (iteration > 2) {
      d.resize(n);
      for (std::size_t i = 0; i < n; ++i) d[i] = a2_[i] - a[i];
      const double d22 = Dot(d, d);
      const double d12 = Dot(c, d);
      const double rd2m = Dot(a, d);
      const double aa2 = d11 * d22;
      const double det = aa2 - d12 * d12;
      double t1, t2;
      if (std::abs(det / aa2) < kDetTolerance) {
        t1 = -rd1m / d11;
        t2 = 0.0;
      } else {
        t1 = (-rd1m * d22 + rd2m * d12) / det;
        t2 = (rd1m * d12 - rd2m * d11) / det;
      }
      for (std::size_t i = 0; i < n; ++i) {
        next[i] = (1.0 - t1 - t2) * b[i] + (t1 * b1_[i] + t2 * b2_[i]) +
                  beta_ * (a[i] + t1 * c[i] + t2 * d[i]);
      }
    } else {
      const double t1 = -rd1m / d11;
      for (std::size_t i = 0; i < n; ++i) {
        next[i] = (1.0 - t1) * b[i] + t1 * b1_[i] + beta_ * (a[i] + t1 * c[i]);
      }
    }
    b2_ = b1_;
    a2_ = a1_;
    a1_ = a;
    b1_ = b;
    b = next;
    return res;
  }

 private:
  static constexpr double kDetTolerance = 1.0e-9;
  double beta_;
  double threshold_;
  std::vector<double> a1_, b1_, a2_, b2_;
};

}  // namespace

std::unique_ptr<Mixer> MakeAndersonMixer(double beta, double threshold) {
  return std::make_unique<AndersonMixer>(beta, threshold);
}

}  // namespace atomx::scf
