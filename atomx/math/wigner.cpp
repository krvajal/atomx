#include "atomx/math/wigner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace atomx::math {

namespace {

double Factorial(int n) { return std::tgamma(static_cast<double>(n) + 1.0); }

bool IsEven(int n) { return n % 2 == 0; }

}  // namespace

double ClebschGordan(int two_j1, int two_m1, int two_j2, int two_m2,
                     int two_J, int two_M) {
  if (two_m1 + two_m2 != two_M) return 0.0;
  if (std::abs(two_m1) > two_j1 || std::abs(two_m2) > two_j2 ||
      std::abs(two_M) > two_J) {
    return 0.0;
  }
  if (!IsEven(two_j1 + two_m1) || !IsEven(two_j2 + two_m2) ||
      !IsEven(two_J + two_M)) {
    return 0.0;
  }
  // Triangle condition, with all combinations below integer.
  const int a = two_J + two_j1 - two_j2;
  const int b = two_J - two_j1 + two_j2;
  const int c = two_j1 + two_j2 - two_J;
  if (a < 0 || b < 0 || c < 0 || !IsEven(a) || !IsEven(b) || !IsEven(c)) {
    return 0.0;
  }
  const int A = a / 2, B = b / 2, C = c / 2;
  const int j1pm1 = (two_j1 + two_m1) / 2, j1mm1 = (two_j1 - two_m1) / 2;
  const int j2pm2 = (two_j2 + two_m2) / 2, j2mm2 = (two_j2 - two_m2) / 2;
  const int JpM = (two_J + two_M) / 2, JmM = (two_J - two_M) / 2;
  const int jsum = (two_j1 + two_j2 + two_J) / 2;

  const double pref =
      std::sqrt((two_J + 1.0) * Factorial(A) * Factorial(B) * Factorial(C) /
                Factorial(jsum + 1)) *
      std::sqrt(Factorial(JpM) * Factorial(JmM) * Factorial(j1mm1) *
                Factorial(j1pm1) * Factorial(j2mm2) * Factorial(j2pm2));

  // Denominators: k, C-k, j1-m1-k, j2+m2-k, J-j2+m1+k, J-j1-m2+k.
  const int d5 = (two_J - two_j2 + two_m1) / 2;
  const int d6 = (two_J - two_j1 - two_m2) / 2;
  const int kmin = std::max({0, -d5, -d6});
  const int kmax = std::min({C, j1mm1, j2pm2});
  double sum = 0.0;
  for (int k = kmin; k <= kmax; ++k) {
    const double sign = IsEven(k) ? 1.0 : -1.0;
    sum += sign / (Factorial(k) * Factorial(C - k) * Factorial(j1mm1 - k) *
                   Factorial(j2pm2 - k) * Factorial(d5 + k) *
                   Factorial(d6 + k));
  }
  return pref * sum;
}

double Wigner3j(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2,
                int two_m3) {
  if (two_m1 + two_m2 + two_m3 != 0) return 0.0;
  const int phase2 = two_j1 - two_j2 - two_m3;  // 2 (j1 - j2 - m3)
  const double sign = IsEven(phase2 / 2) ? 1.0 : -1.0;
  return sign / std::sqrt(two_j3 + 1.0) *
         ClebschGordan(two_j1, two_m1, two_j2, two_m2, two_j3, -two_m3);
}

double LegendreTripleIntegral(int l1, int l2, int l3) {
  const double w = Wigner3j(2 * l1, 2 * l2, 2 * l3, 0, 0, 0);
  return 2.0 * w * w;
}

}  // namespace atomx::math
