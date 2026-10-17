#include "atomx/atom/nucleus.hpp"

#include <cmath>
#include <cstdio>

#include "atomx/atom/periodic_table.hpp"
#include "atomx/core/constants.hpp"

namespace atomx::atom {

using constants::kE2;

Field PointNucleus::Potential(const grid::RadialGrid& grid) const {
  Field v(grid.size());
  for (std::size_t i = 0; i < v.size(); ++i) v[i] = -charge() * kE2 / grid.r(i);
  return v;
}

UniformSphereNucleus::UniformSphereNucleus(double z, double radius)
    : Nucleus(z), radius_(radius) {
  if (radius <= 0.0) throw Error("nuclear radius must be positive");
}

double UniformSphereNucleus::DefaultRadius(int z) {
  const double fermi =
      1.0793 * std::pow(AtomicWeight(z), 1.0 / 3.0) + 0.73587;  // fm
  return fermi * 1.0e-5 / constants::kBohrRadiusAngstrom;
}

std::string UniformSphereNucleus::Describe() const {
  char buf[96];
  std::snprintf(buf, sizeof buf, "uniformly charged sphere, R = %.5e a.u.",
                radius_);
  return buf;
}

Field UniformSphereNucleus::Potential(const grid::RadialGrid& grid) const {
  Field v(grid.size());
  const double R = radius_;
  for (std::size_t i = 0; i < v.size(); ++i) {
    const double r = grid.r(i);
    v[i] = (r < R) ? -charge() * kE2 / (2.0 * R) * (3.0 - (r / R) * (r / R))
                   : -charge() * kE2 / r;
  }
  return v;
}

}  // namespace atomx::atom
