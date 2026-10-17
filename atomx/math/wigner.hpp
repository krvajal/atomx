#pragma once

namespace atomx::math {

// Angular momenta are passed doubled (2j, 2m) so that half-integers are exact.

// Clebsch-Gordan coefficient <j1 m1, j2 m2 | J M> (Racah formula).
double ClebschGordan(int two_j1, int two_m1, int two_j2, int two_m2,
                     int two_J, int two_M);

// Wigner 3j symbol (j1 j2 j3; m1 m2 m3).
double Wigner3j(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2,
                int two_m3);

// \int_{-1}^{1} P_l1(x) P_l2(x) P_l3(x) dx = 2 (l1 l2 l3; 0 0 0)^2,
// the angular factor of the non-relativistic exchange (ld1.x's sl3).
double LegendreTripleIntegral(int l1, int l2, int l3);

}  // namespace atomx::math
