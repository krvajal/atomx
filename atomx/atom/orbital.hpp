#pragma once

#include <string>

namespace atomx::atom {

// One radial orbital (a shell nl, or nlj in the Dirac case) and its
// occupation. For spin-polarised calculations each spin has its own orbitals.
struct Orbital {
  int n = 1;
  int l = 0;
  double j = 0.0;     // total angular momentum; 0 unless relativity is Dirac
  int spin = 0;       // 0 (up / unpolarised) or 1 (down)
  double occupation = 0.0;

  // "1S", "2P", ... as printed by ld1.x
  std::string Label() const;

  bool IsOccupied() const { return occupation > 0.0; }
  // Orbitals with negative occupation are excluded from the calculation.
  bool IsActive() const { return occupation >= 0.0; }

  // Dirac quantum number kappa: -(l+1) for j = l+1/2, l for j = l-1/2.
  int Kappa() const;
};

char AngularMomentumLetter(int l);
int AngularMomentumFromLetter(char c);  // -1 if not s, p, d, f

}  // namespace atomx::atom
