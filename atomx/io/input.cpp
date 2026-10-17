#include "atomx/io/input.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>

#include "atomx/atom/periodic_table.hpp"

namespace atomx::io {

namespace {

std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

}  // namespace

CalculationInput CalculationInput::FromNamelist(const Namelist& nl) {
  CalculationInput in;
  in.title = nl.GetString("title", "");

  const double zed = nl.GetDouble("zed", 0.0);
  const std::string symbol = nl.GetString("atom", "");
  if (zed > 0.0) {
    if (std::abs(zed - std::round(zed)) > 1e-8) {
      throw Error("integer zed is required");
    }
    in.z = static_cast<int>(std::lround(zed));
    if (!symbol.empty() && atom::AtomicNumber(symbol) != in.z) {
      throw Error("zed and atom are inconsistent");
    }
  } else if (!symbol.empty()) {
    in.z = atom::AtomicNumber(symbol);
  } else {
    throw Error("either zed or atom must be given");
  }
  if (in.z < 1 || in.z > atom::kMaxAtomicNumber) {
    throw Error("zed out of range");
  }

  in.configuration = nl.GetString("config", "default");
  if (in.configuration.empty()) {
    throw Error("config is required (reading orbitals from cards is not supported)");
  }

  // ld1.x: rel=5 (default) means 0 below Z=19, scalar-relativistic above.
  int rel = nl.GetInt("rel", 5);
  if (rel == 5) rel = in.z < 19 ? 0 : 1;
  switch (rel) {
    case 0: in.relativity = Relativity::kNone; break;
    case 2: in.relativity = Relativity::kDirac; break;
    case 1:
      throw Error("scalar-relativistic (rel=1) is not implemented; use rel=0 or rel=2");
    default: throw Error("rel must be 0, 1 or 2");
  }
  in.spin_polarized = nl.GetInt("lsd", 0) == 1;
  const std::string dist = Lower(nl.GetString("rel_dist", "energy"));
  if (dist == "average") {
    in.distribution = atom::SpinOrbitDistribution::kAverage;
  } else if (dist != "energy") {
    throw Error("rel_dist must be 'energy' or 'average'");
  }

  in.dft = nl.GetString("dft", "PZ");
  in.xmin = nl.GetDouble("xmin", in.relativity == Relativity::kNone ? -7.0 : -8.0);
  in.dx = nl.GetDouble("dx", 0.008);
  in.rmax = nl.GetDouble("rmax", 100.0);
  if (in.xmin > -2.0) throw Error("wrong xmin");
  if (in.dx <= 0.0) throw Error("wrong dx");

  in.beta = nl.GetDouble("beta", 0.2);
  in.tr2 = nl.GetDouble("tr2", 1e-14);
  in.nucleus = Lower(nl.GetString("nucleus", "point"));
  in.nuclear_radius = nl.GetDouble("rnucleus", 0.0);
  in.speed_of_light = nl.GetDouble("cau_fact", 0.0);
  in.linear_solver = Lower(nl.GetString("linear_solver", "lapack"));

  if (nl.GetInt("iswitch", 1) != 1) {
    throw Error("only all-electron calculations (iswitch=1) are supported");
  }
  if (nl.GetInt("isic", 0) != 0) throw Error("isic is not supported");
  if (nl.GetInt("latt", 0) != 0) throw Error("latt is not supported");
  // Accepted for compatibility with ld1.x inputs, no effect here.
  (void)nl.GetString("prefix", "");
  (void)nl.GetString("file_charge", "");
  (void)nl.GetString("verbosity", "");
  return in;
}

}  // namespace atomx::io
