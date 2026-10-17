#pragma once

// Physical constants, identical to Quantum ESPRESSO's Modules/constants.f90
// (CODATA 2018), so that results can be compared digit by digit with ld1.x.
// Internally all energies are in Rydberg and lengths in bohr, as in ld1.x.

namespace atomx::constants {

inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kFourPi = 4.0 * kPi;

// e^2 in Rydberg atomic units.
inline constexpr double kE2 = 2.0;

inline constexpr double kBohrRadiusSi = 0.529177210903e-10;   // m
inline constexpr double kHartreeSi = 4.3597447222071e-18;     // J
inline constexpr double kElectronVoltSi = 1.602176634e-19;    // J
// Computed as in QE (BOHR_RADIUS_SI * 100 * 1e8) to get the same last bit.
inline constexpr double kBohrRadiusAngstrom = kBohrRadiusSi * 100.0 * 1.0e8;
inline constexpr double kHartreeToEv = kHartreeSi / kElectronVoltSi;
inline constexpr double kRydbergToEv = kHartreeToEv / 2.0;

// Speed of light in Hartree atomic units (QE's C_AU).
inline constexpr double kSpeedOfLight = 137.035999083817785;

}  // namespace atomx::constants
