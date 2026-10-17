#pragma once

#include <string>

namespace atomx::atom {

inline constexpr int kMaxAtomicNumber = 103;

// Element symbol ("H", "He", ...).
std::string ElementSymbol(int z);
// Atomic number from a symbol; throws if unknown.
int AtomicNumber(const std::string& symbol);
// Standard atomic weight (same table as QE's atom_weight).
double AtomicWeight(int z);
// Experimental ground-state configuration, as in ld1.x's default_conf,
// e.g. "[Ar] 4s2.0 3d10.0".
std::string DefaultConfiguration(int z);

}  // namespace atomx::atom
