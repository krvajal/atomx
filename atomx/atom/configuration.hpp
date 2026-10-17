#pragma once

#include <string>
#include <vector>

#include "atomx/atom/orbital.hpp"
#include "atomx/core/types.hpp"

namespace atomx::atom {

// How electrons of an unsplit nl shell are distributed over j = l -/+ 1/2 in
// a Dirac calculation (ld1.x's rel_dist).
enum class SpinOrbitDistribution {
  kEnergy,   // fill j = l-1/2 first
  kAverage,  // proportionally to 2j+1
};

struct ConfigurationOptions {
  Relativity relativity = Relativity::kNone;
  bool spin_polarized = false;
  SpinOrbitDistribution distribution = SpinOrbitDistribution::kEnergy;
};

// The list of orbitals of an atom, built from a configuration string such as
// "[Ar] 4s2 3d10" with the same rules as ld1.x:
//  - spin-polarised: each shell is split into spin-up and spin-down orbitals
//    following Hund's rule (up filled first);
//  - Dirac: each l > 0 shell is split into j = l-1/2 and j = l+1/2.
class Configuration {
 public:
  static Configuration Parse(const std::string& text,
                             const ConfigurationOptions& options);
  // The experimental ground-state configuration of element z.
  static Configuration Default(int z, const ConfigurationOptions& options);

  const std::vector<Orbital>& orbitals() const { return orbitals_; }
  int nspin() const { return nspin_; }
  double ElectronCount() const;

 private:
  std::vector<Orbital> orbitals_;
  int nspin_ = 1;
};

}  // namespace atomx::atom
