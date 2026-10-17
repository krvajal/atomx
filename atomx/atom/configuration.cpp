#include "atomx/atom/configuration.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

#include "atomx/atom/periodic_table.hpp"

namespace atomx::atom {

namespace {

// Closed shells of the noble-gas cores, in ld1.x's order.
struct CoreShell {
  int n, l;
  double occupation;
};
constexpr CoreShell kCoreShells[] = {
    {1, 0, 2},  {2, 0, 2}, {2, 1, 6},  {3, 0, 2}, {3, 1, 6},
    {4, 0, 2},  {4, 1, 6}, {3, 2, 10}, {5, 0, 2}, {5, 1, 6},
    {4, 2, 10}, {6, 0, 2}, {5, 2, 10}, {4, 3, 14}, {6, 1, 6}};

int CoreShellCount(const std::string& core) {
  std::string c = core;
  std::transform(c.begin(), c.end(), c.begin(),
                 [](unsigned char ch) { return std::tolower(ch); });
  if (c == "he") return 1;
  if (c == "ne") return 3;
  if (c == "ar") return 5;
  if (c == "kr") return 8;
  if (c == "xe") return 11;
  if (c == "hg") return 14;
  if (c == "rn") return 15;
  throw Error("unknown core [" + core + "]");
}

std::vector<std::string> Tokenize(const std::string& text) {
  std::string s = text;
  std::replace(s.begin(), s.end(), ',', ' ');
  std::istringstream in(s);
  std::vector<std::string> tokens;
  for (std::string t; in >> t;) tokens.push_back(t);
  return tokens;
}

Orbital ParseShell(const std::string& token) {
  if (token.size() < 3 || !std::isdigit(static_cast<unsigned char>(token[0]))) {
    throw Error("cannot parse orbital '" + token + "'");
  }
  Orbital o;
  o.n = token[0] - '0';
  o.l = AngularMomentumFromLetter(token[1]);
  if (o.n <= 0 || o.n > 7) throw Error("wrong main quantum number: " + token);
  if (o.l < 0) throw Error("l not found: " + token);
  if (o.l >= o.n) throw Error("main/angular quantum number mismatch: " + token);
  try {
    o.occupation = std::stod(token.substr(2));
  } catch (const std::exception&) {
    throw Error("no occupation in '" + token + "'");
  }
  if (o.occupation > 2.0 * (2 * o.l + 1)) {
    throw Error("wrong occupation: " + token);
  }
  return o;
}

bool SameShell(const Orbital& a, const Orbital& b) {
  return a.n == b.n && a.l == b.l;
}

// el_config: the explicit list, marking repeated shells as the second spin or
// the second j.
std::vector<Orbital> ParseList(const std::string& text,
                               const ConfigurationOptions& opt) {
  auto tokens = Tokenize(text);
  if (tokens.empty()) throw Error("empty configuration");
  std::vector<Orbital> orbitals;
  std::size_t first = 0;
  const std::string& t0 = tokens[0];
  if (t0.size() >= 4 && t0.front() == '[' && t0[3] == ']') {
    const int ncore = CoreShellCount(t0.substr(1, 2));
    for (int i = 0; i < ncore; ++i) {
      Orbital o;
      o.n = kCoreShells[i].n;
      o.l = kCoreShells[i].l;
      o.occupation = kCoreShells[i].occupation;
      orbitals.push_back(o);
    }
    // "[Ar]4s2" without a space
    if (t0.size() > 4) {
      tokens[0] = t0.substr(4);
    } else {
      first = 1;
    }
  }
  const std::size_t ncore = orbitals.size();
  for (std::size_t i = first; i < tokens.size(); ++i) {
    orbitals.push_back(ParseShell(tokens[i]));
  }
  for (std::size_t k = ncore; k < orbitals.size(); ++k) {
    Orbital& o = orbitals[k];
    for (std::size_t i = 0; i < k; ++i) {
      Orbital& p = orbitals[i];
      if (!SameShell(p, o)) continue;
      bool too_many = true;
      if (opt.relativity == Relativity::kDirac) {
        too_many = p.j > 0.0 || o.l == 0;
        p.j = p.l - 0.5;
        o.j = o.l + 0.5;
        if (o.occupation > 2 * o.j + 1 || p.occupation > 2 * p.j + 1) {
          throw Error("occupation wrong for " + o.Label());
        }
      } else if (opt.spin_polarized) {
        too_many = p.spin > 0;
        o.spin = 1;
      }
      if (too_many) throw Error(o.Label() + " found too many times");
    }
  }
  return orbitals;
}

// occ_spin: Hund's rule split into up and down orbitals.
void SplitSpin(std::vector<Orbital>& orbitals) {
  const std::size_t n0 = orbitals.size();
  for (std::size_t n = 0; n < n0; ++n) {
    Orbital& o = orbitals[n];
    const double half = 2 * o.l + 1;
    Orbital partner = o;
    partner.spin = 1 - o.spin;
    if (o.occupation > half) {
      for (std::size_t m = n + 1; m < n0; ++m) {
        if (SameShell(orbitals[m], o)) throw Error("wrong occupations");
      }
      partner.occupation = o.occupation - half;
      orbitals[n].occupation = half;
      orbitals.push_back(partner);
    } else {
      bool unique = true;
      for (std::size_t m = 0; m < n0; ++m) {
        if (m != n && SameShell(orbitals[m], o)) unique = false;
      }
      if (unique) {
        partner.occupation = o.occupation < 0.0 ? o.occupation : 0.0;
        orbitals.push_back(partner);
      }
    }
  }
}

// occ_spinorb: j-split every l > 0 shell given without j.
void SplitSpinOrbit(std::vector<Orbital>& orbitals,
                    SpinOrbitDistribution dist) {
  std::vector<Orbital> out;
  for (Orbital o : orbitals) {
    if (o.l > 0 && std::abs(o.j) < 1e-2) {
      Orbital low = o, high = o;
      low.j = o.l - 0.5;
      high.j = o.l + 0.5;
      if (o.occupation > -1e-3) {
        if (dist == SpinOrbitDistribution::kAverage) {
          high.occupation = o.occupation * (2.0 * (o.l + 1)) / (2.0 * (2 * o.l + 1));
          low.occupation = o.occupation * (2.0 * o.l) / (2.0 * (2 * o.l + 1));
        } else {
          low.occupation = std::min(2.0 * o.l, o.occupation);
          high.occupation =
              std::max(0.0, o.occupation - std::min(2.0 * o.l, o.occupation));
        }
      }
      out.push_back(low);
      out.push_back(high);
    } else {
      if (std::abs(o.j) < 1e-2) o.j = 0.5;
      out.push_back(o);
    }
  }
  orbitals = std::move(out);
}

}  // namespace

Configuration Configuration::Parse(const std::string& text,
                                   const ConfigurationOptions& options) {
  if (options.relativity == Relativity::kDirac && options.spin_polarized) {
    throw Error("spin polarisation and spin-orbit are not allowed together");
  }
  Configuration c;
  c.orbitals_ = ParseList(text, options);
  if (options.spin_polarized) {
    SplitSpin(c.orbitals_);
    c.nspin_ = 2;
  } else if (options.relativity == Relativity::kDirac) {
    SplitSpinOrbit(c.orbitals_, options.distribution);
  }
  return c;
}

Configuration Configuration::Default(int z,
                                     const ConfigurationOptions& options) {
  return Parse(DefaultConfiguration(z), options);
}

double Configuration::ElectronCount() const {
  double n = 0.0;
  for (const auto& o : orbitals_) n += std::max(0.0, o.occupation);
  return n;
}

}  // namespace atomx::atom
