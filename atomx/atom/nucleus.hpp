#pragma once

#include <memory>
#include <string>

#include "atomx/core/types.hpp"
#include "atomx/grid/radial_grid.hpp"

namespace atomx::atom {

// The electron-nucleus interaction.
class Nucleus {
 public:
  virtual ~Nucleus() = default;

  virtual std::string Describe() const = 0;
  double charge() const { return z_; }
  // true if the potential is -Z e^2 / r all the way to the origin; radial
  // solvers use this to pick their series start at r -> 0.
  virtual bool IsPointLike() const = 0;
  // Potential energy of an electron, in Ry.
  virtual Field Potential(const grid::RadialGrid& grid) const = 0;

 protected:
  explicit Nucleus(double z) : z_(z) {}

 private:
  double z_;
};

class PointNucleus final : public Nucleus {
 public:
  explicit PointNucleus(double z) : Nucleus(z) {}
  std::string Describe() const override { return "point nucleus"; }
  bool IsPointLike() const override { return true; }
  Field Potential(const grid::RadialGrid& grid) const override;
};

// Uniformly charged sphere of radius R:
//   V(r) = -Z e^2 / (2R) (3 - r^2/R^2) for r < R, -Z e^2 / r outside.
class UniformSphereNucleus final : public Nucleus {
 public:
  UniformSphereNucleus(double z, double radius);
  // R = 1.0793 A^(1/3) + 0.73587 fm with A the standard atomic weight
  // (Engel et al., Phys. Rev. A 58, 964 (1998)).
  static double DefaultRadius(int z);

  std::string Describe() const override;
  bool IsPointLike() const override { return false; }
  Field Potential(const grid::RadialGrid& grid) const override;
  double radius() const { return radius_; }

 private:
  double radius_;
};

}  // namespace atomx::atom
