#include "atomx/atom/orbital.hpp"

#include <cctype>
#include <cmath>

namespace atomx::atom {

char AngularMomentumLetter(int l) {
  static constexpr char kLetters[] = "SPDFGHI";
  return (l >= 0 && l < 7) ? kLetters[l] : '?';
}

int AngularMomentumFromLetter(char c) {
  switch (std::toupper(static_cast<unsigned char>(c))) {
    case 'S': return 0;
    case 'P': return 1;
    case 'D': return 2;
    case 'F': return 3;
    default: return -1;
  }
}

std::string Orbital::Label() const {
  return std::to_string(n) + AngularMomentumLetter(l);
}

int Orbital::Kappa() const {
  return (std::abs(j - (l + 0.5)) < 1e-6) ? -(l + 1) : l;
}

}  // namespace atomx::atom
