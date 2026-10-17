#include "atomx/atom/periodic_table.hpp"

#include <array>
#include <string_view>

#include "atomx/core/types.hpp"

namespace atomx::atom {

namespace {

// Generated from Quantum ESPRESSO 7.6: atomic/src/default_conf.f90 and
// Modules/atom_weight.f90.
struct Element {
  std::string_view symbol;
  double weight;
  std::string_view configuration;
};

constexpr std::array<Element, kMaxAtomicNumber> kElements = {{
    {"H", 1.00794, "1s1.0"},
    {"He", 4.0026, "1s2.0"},
    {"Li", 6.941, "[He] 2s1.0"},
    {"Be", 9.01218, "[He] 2s2.0"},
    {"B", 10.811, "[He] 2s2.0 2p1.0"},
    {"C", 12.0107, "[He] 2s2.0 2p2.0"},
    {"N", 14.00674, "[He] 2s2.0 2p3.0"},
    {"O", 15.9994, "[He] 2s2.0 2p4.0"},
    {"F", 18.9984, "[He] 2s2.0 2p5.0"},
    {"Ne", 20.1797, "[He] 2s2.0 2p6.0"},
    {"Na", 22.98977, "[Ne] 3s1.0"},
    {"Mg", 24.305, "[Ne] 3s2.0"},
    {"Al", 26.98154, "[Ne] 3s2.0 3p1.0"},
    {"Si", 28.0855, "[Ne] 3s2.0 3p2.0"},
    {"P", 30.97376, "[Ne] 3s2.0 3p3.0"},
    {"S", 32.066, "[Ne] 3s2.0 3p4.0"},
    {"Cl", 35.4527, "[Ne] 3s2.0 3p5.0"},
    {"Ar", 39.948, "[Ne] 3s2.0 3p6.0"},
    {"K", 39.0983, "[Ar] 4s1.0"},
    {"Ca", 40.078, "[Ar] 4s2.0"},
    {"Sc", 44.95591, "[Ar] 4s2.0 3d1.0"},
    {"Ti", 47.867, "[Ar] 4s2.0 3d2.0"},
    {"V", 50.9415, "[Ar] 4s2.0 3d3.0"},
    {"Cr", 51.9961, "[Ar] 4s1.0 3d5.0"},
    {"Mn", 54.93805, "[Ar] 4s2.0 3d5.0"},
    {"Fe", 55.845, "[Ar] 4s2.0 3d6.0"},
    {"Co", 58.9332, "[Ar] 4s2.0 3d7.0"},
    {"Ni", 58.6934, "[Ar] 4s2.0 3d8.0"},
    {"Cu", 63.546, "[Ar] 4s1.0 3d10.0"},
    {"Zn", 65.39, "[Ar] 4s2.0 3d10.0"},
    {"Ga", 69.723, "[Ar] 4s2.0 3d10.0 4p1.0"},
    {"Ge", 72.61, "[Ar] 4s2.0 3d10.0 4p2.0"},
    {"As", 74.9216, "[Ar] 4s2.0 3d10.0 4p3.0"},
    {"Se", 78.96, "[Ar] 4s2.0 3d10.0 4p4.0"},
    {"Br", 79.904, "[Ar] 4s2.0 3d10.0 4p5.0"},
    {"Kr", 83.8, "[Ar] 4s2.0 3d10.0 4p6.0"},
    {"Rb", 85.4678, "[Kr] 5s1.0"},
    {"Sr", 87.62, "[Kr] 5s2.0"},
    {"Y", 88.90585, "[Kr] 5s2.0 4d1.0"},
    {"Zr", 91.224, "[Kr] 5s2.0 4d2.0"},
    {"Nb", 92.90638, "[Kr] 5s1.0 4d4.0"},
    {"Mo", 95.94, "[Kr] 5s1.0 4d5.0"},
    {"Tc", 98.0, "[Kr] 5s2.0 4d5.0"},
    {"Ru", 101.07, "[Kr] 5s1.0 4d7.0"},
    {"Rh", 102.9055, "[Kr] 5s1.0 4d8.0"},
    {"Pd", 106.42, "[Kr] 5s0.0 4d10.0"},
    {"Ag", 107.8682, "[Kr] 5s1.0 4d10.0"},
    {"Cd", 112.411, "[Kr] 5s2.0 4d10.0"},
    {"In", 114.818, "[Kr] 5s2.0 4d10.0 5p1.0"},
    {"Sn", 118.71, "[Kr] 5s2.0 4d10.0 5p2.0"},
    {"Sb", 121.76, "[Kr] 5s2.0 4d10.0 5p3.0"},
    {"Te", 127.6, "[Kr] 5s2.0 4d10.0 5p4.0"},
    {"I", 126.90447, "[Kr] 5s2.0 4d10.0 5p5.0"},
    {"Xe", 131.29, "[Kr] 5s2.0 4d10.0 5p6.0"},
    {"Cs", 132.90545, "[Xe] 6s1.0"},
    {"Ba", 137.327, "[Xe] 6s2.0"},
    {"La", 138.9055, "[Xe] 6s2.0 5d1.0 4f0.0"},
    {"Ce", 140.116, "[Xe] 6s2.0 5d1.0 4f1.0"},
    {"Pr", 140.90765, "[Xe] 6s2.0 5d0.0 4f3.0"},
    {"Nd", 144.24, "[Xe] 6s2.0 5d0.0 4f4.0"},
    {"Pm", 145.0, "[Xe] 6s2.0 5d0.0 4f5.0"},
    {"Sm", 150.36, "[Xe] 6s2.0 5d0.0 4f6.0"},
    {"Eu", 151.964, "[Xe] 6s2.0 5d0.0 4f7.0"},
    {"Gd", 157.25, "[Xe] 6s2.0 5d1.0 4f7.0"},
    {"Tb", 158.92534, "[Xe] 6s2.0 5d0.0 4f9.0"},
    {"Dy", 162.5, "[Xe] 6s2.0 5d0.0 4f10.0"},
    {"Ho", 164.93032, "[Xe] 6s2.0 5d0.0 4f11.0"},
    {"Er", 167.26, "[Xe] 6s2.0 5d0.0 4f12.0"},
    {"Tm", 168.93421, "[Xe] 6s2.0 5d0.0 4f13.0"},
    {"Yb", 173.04, "[Xe] 6s2.0 5d0.0 4f14.0"},
    {"Lu", 174.967, "[Xe] 6s2.0 5d1.0 4f14.0"},
    {"Hf", 178.49, "[Xe] 6s2.0 5d2.0 4f14.0"},
    {"Ta", 180.9479, "[Xe] 6s2.0 5d3.0 4f14.0"},
    {"W", 183.84, "[Xe] 6s2.0 5d4.0 4f14.0"},
    {"Re", 186.207, "[Xe] 6s2.0 5d5.0 4f14.0"},
    {"Os", 190.23, "[Xe] 6s2.0 5d6.0 4f14.0"},
    {"Ir", 192.217, "[Xe] 6s2.0 5d7.0 4f14.0"},
    {"Pt", 195.078, "[Xe] 6s1.0 5d9.0 4f14.0"},
    {"Au", 196.96655, "[Xe] 6s1.0 5d10.0 4f14.0"},
    {"Hg", 200.59, "[Xe] 6s2.0 5d10.0 4f14.0"},
    {"Tl", 204.3833, "[Xe] 6s2.0 5d10.0 4f14.0 6p1.0"},
    {"Pb", 207.2, "[Xe] 6s2.0 5d10.0 4f14.0 6p2.0"},
    {"Bi", 208.98038, "[Xe] 6s2.0 5d10.0 4f14.0 6p3.0"},
    {"Po", 209.0, "[Xe] 6s2.0 5d10.0 4f14.0 6p4.0"},
    {"At", 210.0, "[Xe] 6s2.0 5d10.0 4f14.0 6p5.0"},
    {"Rn", 222.0, "[Xe] 6s2.0 5d10.0 4f14.0 6p6.0"},
    {"Fr", 223.0, "[Rn] 7s1.0"},
    {"Ra", 226.0, "[Rn] 7s2.0"},
    {"Ac", 227.0, "[Rn] 7s2.0 6d1.0 5f0.0"},
    {"Th", 232.0381, "[Rn] 7s2.0 6d2.0 5f0.0"},
    {"Pa", 231.03588, "[Rn] 7s2.0 6d1.0 5f2.0"},
    {"U", 238.0289, "[Rn] 7s2.0 6d1.0 5f3.0"},
    {"Np", 237.0, "[Rn] 7s2.0 6d1.0 5f4.0"},
    {"Pu", 244.0, "[Rn] 7s2.0 6d0.0 5f6.0"},
    {"Am", 243.0, "[Rn] 7s2.0 6d0.0 5f7.0"},
    {"Cm", 247.0, "[Rn] 7s2.0 6d1.0 5f7.0"},
    {"Bk", 247.0, "[Rn] 7s2.0 6d0.0 5f9.0"},
    {"Cf", 251.0, "[Rn] 7s2.0 6d0.0 5f10.0"},
    {"Es", 252.0, "[Rn] 7s2.0 6d0.0 5f11.0"},
    {"Fm", 257.0, "[Rn] 7s2.0 6d0.0 5f12.0"},
    {"Md", 258.0, "[Rn] 7s2.0 6d0.0 5f13.0"},
    {"No", 259.0, "[Rn] 7s2.0 6d0.0 5f14.0"},
    {"Lr", 262.0, "[Rn] 7s2.0 7p1.0 5f14.0"},
}};

const Element& Lookup(int z) {
  if (z < 1 || z > kMaxAtomicNumber) {
    throw Error("atomic number out of range: " + std::to_string(z));
  }
  return kElements[z - 1];
}

}  // namespace

std::string ElementSymbol(int z) { return std::string(Lookup(z).symbol); }

int AtomicNumber(const std::string& symbol) {
  for (int z = 1; z <= kMaxAtomicNumber; ++z) {
    if (kElements[z - 1].symbol == symbol) return z;
  }
  throw Error("unknown element symbol: " + symbol);
}

double AtomicWeight(int z) { return Lookup(z).weight; }

std::string DefaultConfiguration(int z) {
  return std::string(Lookup(z).configuration);
}

}  // namespace atomx::atom
