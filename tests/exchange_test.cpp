// Exchange operators and local exchange-correlation functionals against
// analytic results, on fixed (non-self-consistent) orbitals.

#include <cmath>

#include <gtest/gtest.h>

#include "atomx/atom/atom_state.hpp"
#include "atomx/atom/configuration.hpp"
#include "atomx/core/constants.hpp"
#include "atomx/grid/poisson.hpp"
#include "atomx/grid/radial_grid.hpp"
#include "atomx/xc/exchange_operator.hpp"
#include "atomx/xc/lda.hpp"
#include "tests/test_util.hpp"

namespace atomx {
namespace {

using testing::Hydrogen1s;
using testing::Hydrogen2p;
using testing::Hydrogen2s;
using testing::SetOrbital;

double ExchangeEnergy(const xc::ExchangeOperator& op,
                      const atom::AtomState& state) {
  return xc::ExchangeOperator::Energy(state, op.ApplyAll(state));
}

TEST(ExchangeOperator, HydrogenicHeliumIsMinusJ) {
  // Two electrons in one hydrogenic 1s orbital: E_x = -J = -5Z/8 Ha.
  const double z = 2.0;
  const auto grid = grid::RadialGrid::Logarithmic(-9.0, 0.004, 60.0, z);
  const grid::RadialPoisson poisson(grid);
  const auto nonrel = xc::MakeNonRelativisticExchange(poisson);
  const auto dirac = xc::MakeDiracExchange(poisson);

  atom::AtomState s0(grid, atom::Configuration::Parse("1s2", {}),
                     Relativity::kNone);
  SetOrbital(s0, 0, [&](double r) { return Hydrogen1s(z, r); });
  EXPECT_NEAR(0.5 * ExchangeEnergy(*nonrel, s0), -5.0 * z / 8.0, 1e-9);

  // The Dirac operator with a zero small component must agree.
  atom::ConfigurationOptions opt;
  opt.relativity = Relativity::kDirac;
  atom::AtomState s2(grid, atom::Configuration::Parse("1s2", opt),
                     Relativity::kDirac);
  SetOrbital(s2, 0, [&](double r) { return Hydrogen1s(z, r); });
  EXPECT_NEAR(0.5 * ExchangeEnergy(*dirac, s2), -5.0 * z / 8.0, 1e-9);
}

TEST(ExchangeOperator, DiracRecouplingOfClosedShellsMatchesNonRelativistic) {
  // With identical radial functions and F = 0, the j-coupled exchange of
  // closed nl shells (p1/2 + p3/2) must equal the l-coupled one. This checks
  // every angular factor, including the cross terms between j shells.
  const double z = 10.0;
  const auto grid = grid::RadialGrid::Logarithmic(-8.0, 0.006, 60.0, z);
  const grid::RadialPoisson poisson(grid);

  atom::AtomState nr(grid, atom::Configuration::Parse("1s2 2s2 2p6 3d10", {}),
                     Relativity::kNone);
  auto shape = [&](int n, int l) {
    return [=](double r) {
      if (l == 0) return n == 1 ? Hydrogen1s(z, r) : Hydrogen2s(z, r);
      if (l == 1) return Hydrogen2p(z, r);
      return r * r * r * std::exp(-z * r / 3.0);  // 3d
    };
  };
  for (std::size_t k = 0; k < nr.size(); ++k) {
    SetOrbital(nr, k, shape(nr.orbitals()[k].n, nr.orbitals()[k].l));
  }

  atom::ConfigurationOptions opt;
  opt.relativity = Relativity::kDirac;
  atom::AtomState dr(grid, atom::Configuration::Parse("1s2 2s2 2p6 3d10", opt),
                     Relativity::kDirac);
  ASSERT_EQ(dr.size(), 6u);  // 1s, 2s, 2p1/2, 2p3/2, 3d3/2, 3d5/2
  for (std::size_t k = 0; k < dr.size(); ++k) {
    SetOrbital(dr, k, shape(dr.orbitals()[k].n, dr.orbitals()[k].l));
  }

  const double e_nr =
      ExchangeEnergy(*xc::MakeNonRelativisticExchange(poisson), nr);
  const double e_dr = ExchangeEnergy(*xc::MakeDiracExchange(poisson), dr);
  EXPECT_NEAR(e_dr, e_nr, 1e-10 * std::abs(e_nr));
}

TEST(ExchangeOperator, OpenShellAverageRemovesSelfInteraction) {
  // One electron in a p shell: exchange must cancel its Hartree energy,
  // E_x = -E_H. This needs the open-shell average term of ld1.x's dvex.
  const double z = 3.0;
  const auto grid = grid::RadialGrid::Logarithmic(-8.0, 0.006, 60.0, z);
  const grid::RadialPoisson poisson(grid);
  // spin-polarised, so that the electron is not split over two spins
  atom::ConfigurationOptions opt;
  opt.spin_polarized = true;
  atom::AtomState s(grid, atom::Configuration::Parse("2p1", opt),
                    Relativity::kNone);
  SetOrbital(s, 0, [&](double r) { return Hydrogen2p(z, r); });
  const double ex = ExchangeEnergy(*xc::MakeNonRelativisticExchange(poisson), s);
  const Field rho = s.Density()[0];  // all in spin up
  const Field vh = poisson.Solve(0, 2, rho);
  Field f(grid.size());
  for (std::size_t i = 0; i < f.size(); ++i) f[i] = rho[i] * vh[i];
  const double eh = 0.5 * constants::kE2 * grid.Integrate(f, 4);
  EXPECT_NEAR(ex, -eh, 1e-9);
}

TEST(LocalDensity, FullyPolarisedSlaterExchangeScalesAsCubeRootOfTwo) {
  // For the same density, a fully spin-polarised electron gas has
  // 2^(1/3) times the exchange energy of an unpolarised one.
  const double z = 1.0;
  const auto grid = grid::RadialGrid::Logarithmic(-7.0, 0.008, 60.0, z);

  atom::AtomState unpol(grid, atom::Configuration::Parse("1s1", {}),
                        Relativity::kNone);
  SetOrbital(unpol, 0, [&](double r) { return Hydrogen1s(z, r); });
  atom::ConfigurationOptions opt;
  opt.spin_polarized = true;
  atom::AtomState pol(grid, atom::Configuration::Parse("1s1", opt),
                      Relativity::kNone);
  ASSERT_EQ(pol.nspin(), 2);
  SetOrbital(pol, 0, [&](double r) { return Hydrogen1s(z, r); });

  xc::SlaterExchange slater;
  const double e_u = slater.Evaluate(unpol, unpol.Density()).energy;
  const double e_p = slater.Evaluate(pol, pol.Density()).energy;
  EXPECT_NEAR(e_p / e_u, std::cbrt(2.0), 1e-12);
  // Hydrogen 1s, unpolarised LDA exchange: -(3/4)(3/pi)^(1/3) \int n^(4/3)
  // = -0.2127 Ha for the exact density (Ry here).
  EXPECT_NEAR(0.5 * e_u, -0.212742, 1e-5);
}

TEST(LocalDensity, PolarisedCorrelationIsWeakerAndSpinPotentialsSplit) {
  const double z = 1.0;
  const auto grid = grid::RadialGrid::Logarithmic(-7.0, 0.008, 60.0, z);
  atom::ConfigurationOptions opt;
  opt.spin_polarized = true;
  atom::AtomState pol(grid, atom::Configuration::Parse("1s1", opt),
                      Relativity::kNone);
  SetOrbital(pol, 0, [&](double r) { return Hydrogen1s(z, r); });
  atom::AtomState unpol(grid, atom::Configuration::Parse("1s1", {}),
                        Relativity::kNone);
  SetOrbital(unpol, 0, [&](double r) { return Hydrogen1s(z, r); });

  xc::PerdewZungerCorrelation pz;
  const auto cp = pz.Evaluate(pol, pol.Density());
  const auto cu = pz.Evaluate(unpol, unpol.Density());
  EXPECT_LT(cu.energy, cp.energy);  // polarised gas: less correlation
  EXPECT_LT(cp.energy, 0.0);
  // spin-up and spin-down potentials differ where the density is polarised
  const std::size_t i = grid.size() / 2;
  EXPECT_GT(std::abs(cp.potential[0][i] - cp.potential[1][i]), 1e-4);
}

}  // namespace
}  // namespace atomx
