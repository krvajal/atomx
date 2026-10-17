#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace atomx {

// A function sampled on the radial grid.
using Field = std::vector<double>;

// One field per spin channel. Only the first `nspin` entries are used.
using SpinField = std::array<Field, 2>;

// Radial two-component spinor (large, small). For non-relativistic and
// scalar-relativistic orbitals the small component is identically zero.
struct Spinor {
  Field large;
  Field small;
};

// How the kinetic energy is treated.
enum class Relativity {
  kNone,    // Schroedinger equation
  kScalar,  // scalar-relativistic (not implemented yet)
  kDirac,   // Dirac equation, j-split orbitals
};

class Error : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

}  // namespace atomx
