// Unit tests of the numerical building blocks against analytic results.

#include <cmath>

#include <gtest/gtest.h>

#include "atomx/atom/configuration.hpp"
#include "atomx/atom/nucleus.hpp"
#include "atomx/core/constants.hpp"
#include "atomx/grid/poisson.hpp"
#include "atomx/grid/radial_grid.hpp"
#include "atomx/math/wigner.hpp"
#include "atomx/solvers/radial_solver.hpp"

namespace atomx {
namespace {

grid::RadialGrid MakeGrid(double z = 1.0) {
  return grid::RadialGrid::Logarithmic(-8.0, 0.005, 80.0, z);
}

TEST(RadialGrid, IntegratesPolynomialTimesExponential) {
  const auto g = MakeGrid();
  Field f(g.size());
  for (std::size_t i = 0; i < f.size(); ++i) {
    f[i] = g.r(i) * g.r(i) * std::exp(-g.r(i));
  }
  EXPECT_NEAR(g.Integrate(f, 2), 2.0, 1e-10);  // Gamma(3)
}

TEST(RadialPoisson, HydrogenGroundStateDensity) {
  const auto g = MakeGrid();
  const grid::RadialPoisson poisson(g);
  Field f(g.size());
  for (std::size_t i = 0; i < f.size(); ++i) {
    const double r = g.r(i);
    f[i] = 4.0 * r * r * std::exp(-2.0 * r);  // 4 pi r^2 |psi_1s|^2
  }
  const Field v = poisson.Solve(0, 2, f);
  for (std::size_t i = 0; i < g.size(); i += 97) {
    const double r = g.r(i);
    const double exact = (1.0 - (1.0 + r) * std::exp(-2.0 * r)) / r;
    EXPECT_NEAR(v[i], exact, 1e-8 * std::max(1.0, exact)) << "r=" << r;
  }
}

TEST(Wigner, KnownCoefficients) {
  // <1/2 1/2, 1/2 -1/2 | 0 0> = 1/sqrt(2)
  EXPECT_NEAR(math::ClebschGordan(1, 1, 1, -1, 0, 0), 1.0 / std::sqrt(2.0), 1e-14);
  // <1 0, 1 0 | 2 0> = sqrt(2/3)
  EXPECT_NEAR(math::ClebschGordan(2, 0, 2, 0, 4, 0), std::sqrt(2.0 / 3.0), 1e-14);
  // \int P_l P_l P_0 = 2/(2l+1)
  for (int l = 0; l <= 6; ++l) {
    EXPECT_NEAR(math::LegendreTripleIntegral(l, l, 0), 2.0 / (2 * l + 1), 1e-14);
  }
  // parity: odd l1+l2+l3 vanishes
  EXPECT_NEAR(math::LegendreTripleIntegral(1, 1, 1), 0.0, 1e-15);
}

TEST(Configuration, CoreAndValence) {
  const auto c = atom::Configuration::Parse("[Ar] 4s2 3d10", {});
  EXPECT_EQ(c.orbitals().size(), 7u);
  EXPECT_DOUBLE_EQ(c.ElectronCount(), 30.0);
  EXPECT_EQ(c.orbitals().back().Label(), "3D");
}

TEST(Configuration, HundsRuleSpinSplit) {
  atom::ConfigurationOptions opt;
  opt.spin_polarized = true;
  const auto c = atom::Configuration::Parse("1s2 2s2 2p3", opt);
  ASSERT_EQ(c.nspin(), 2);
  double up = 0.0, down = 0.0;
  for (const auto& o : c.orbitals()) (o.spin == 0 ? up : down) += o.occupation;
  EXPECT_DOUBLE_EQ(up, 5.0);  // N: 1s1 2s1 2p3 up
  EXPECT_DOUBLE_EQ(down, 2.0);
}

TEST(Configuration, SpinOrbitSplitFillsLowerJFirst) {
  atom::ConfigurationOptions opt;
  opt.relativity = Relativity::kDirac;
  const auto c = atom::Configuration::Parse("1s2 2s2 2p4", opt);
  ASSERT_EQ(c.orbitals().size(), 4u);
  EXPECT_DOUBLE_EQ(c.orbitals()[2].j, 0.5);
  EXPECT_DOUBLE_EQ(c.orbitals()[2].occupation, 2.0);
  EXPECT_DOUBLE_EQ(c.orbitals()[3].j, 1.5);
  EXPECT_DOUBLE_EQ(c.orbitals()[3].occupation, 2.0);
}

TEST(SchrodingerSolver, HydrogenSpectrum) {
  const auto g = MakeGrid();
  const atom::PointNucleus nucleus(1.0);
  const auto solver = solvers::MakeSchrodingerSolver(g, nucleus);
  const Field v = nucleus.Potential(g);
  for (int n = 1; n <= 3; ++n) {
    for (int l = 0; l < n; ++l) {
      atom::Orbital o;
      o.n = n;
      o.l = l;
      const auto sol = solver->Solve(o, v, -1.0 / (n * n) * 0.9);
      ASSERT_TRUE(sol.converged);
      EXPECT_NEAR(sol.energy, -1.0 / (n * n), 1e-8) << o.Label();  // Ry
    }
  }
}

TEST(DiracSolver, HydrogenLikeMercury) {
  const double z = 80.0;
  const auto g = MakeGrid(z);
  const atom::PointNucleus nucleus(z);
  const double c = constants::kSpeedOfLight;
  const auto solver = solvers::MakeDiracSolver(g, nucleus, c);
  atom::Orbital o;  // 1s1/2
  o.j = 0.5;
  const auto sol = solver->Solve(o, nucleus.Potential(g), -z * z);
  ASSERT_TRUE(sol.converged);
  // Exact Dirac energy without rest mass, in Ry.
  const double za = z / c;
  const double exact = 2.0 * c * c * (std::sqrt(1.0 - za * za) - 1.0);
  EXPECT_NEAR(sol.energy, exact, 1e-6 * std::abs(exact));
}

// Exact Dirac energy (Ry, without rest mass) of a hydrogen-like ion.
double DiracEnergy(double z, double c, int n, int kappa) {
  const double za = z / c;
  const double gamma = std::sqrt(kappa * kappa - za * za);
  const double d = n - std::abs(kappa) + gamma;
  return 2.0 * c * c * (1.0 / std::sqrt(1.0 + za * za / (d * d)) - 1.0);
}

TEST(DiracSolver, HydrogenLikeFineStructure) {
  // 2s1/2 and 2p1/2 are degenerate; 2p3/2 lies above. Covers both signs of
  // kappa in the point-nucleus start at the origin.
  const double z = 80.0;
  const auto g = MakeGrid(z);
  const atom::PointNucleus nucleus(z);
  const double c = constants::kSpeedOfLight;
  const auto solver = solvers::MakeDiracSolver(g, nucleus, c);
  const Field v = nucleus.Potential(g);
  struct State { int n, l; double j; int kappa; };
  for (const State st : {State{2, 0, 0.5, -1}, State{2, 1, 0.5, 1},
                         State{2, 1, 1.5, -2}, State{3, 2, 2.5, -3},
                         State{3, 2, 1.5, 2}}) {
    atom::Orbital o;
    o.n = st.n;
    o.l = st.l;
    o.j = st.j;
    ASSERT_EQ(o.Kappa(), st.kappa);
    const auto sol = solver->Solve(o, v, -z * z / (st.n * st.n));
    ASSERT_TRUE(sol.converged) << o.Label() << " j=" << o.j;
    const double exact = DiracEnergy(z, c, st.n, st.kappa);
    EXPECT_NEAR(sol.energy, exact, 1e-7 * std::abs(exact))
        << o.Label() << " j=" << o.j;
  }
}

TEST(DiracSolver, FiniteNucleusRaisesDeepLevels) {
  // A finite nucleus weakens the attraction inside R: 1s moves up, a 2p3/2
  // orbital (which vanishes at the origin) is essentially unaffected.
  const double z = 92.0;
  const auto g = MakeGrid(z);
  const atom::PointNucleus point(z);
  const atom::UniformSphereNucleus sphere(
      z, atom::UniformSphereNucleus::DefaultRadius(92));
  const double c = constants::kSpeedOfLight;
  atom::Orbital s1;
  s1.j = 0.5;
  atom::Orbital p32;
  p32.n = 2;
  p32.l = 1;
  p32.j = 1.5;
  const auto ep = solvers::MakeDiracSolver(g, point, c);
  const auto es = solvers::MakeDiracSolver(g, sphere, c);
  const double s_point = ep->Solve(s1, point.Potential(g), -z * z).energy;
  const double s_sphere = es->Solve(s1, sphere.Potential(g), -z * z).energy;
  // U^91+ 1s finite-size shift: about 198 eV (Ry here)
  EXPECT_NEAR((s_sphere - s_point) * constants::kRydbergToEv, 198.0, 10.0);
  const double p_point = ep->Solve(p32, point.Potential(g), -z * z / 4).energy;
  const double p_sphere = es->Solve(p32, sphere.Potential(g), -z * z / 4).energy;
  EXPECT_NEAR(p_sphere, p_point, 1e-6 * std::abs(p_point));
}

TEST(SchrodingerSolver, RecoversFromPoorEnergyGuess) {
  const auto g = MakeGrid();
  const atom::PointNucleus nucleus(1.0);
  const auto solver = solvers::MakeSchrodingerSolver(g, nucleus);
  atom::Orbital o;
  o.n = 3;
  // far too deep: too few nodes at first
  auto sol = solver->Solve(o, nucleus.Potential(g), -50.0);
  ASSERT_TRUE(sol.converged);
  EXPECT_NEAR(sol.energy, -1.0 / 9.0, 1e-8);
  // far too shallow: too many nodes at first
  sol = solver->Solve(o, nucleus.Potential(g), -1e-4);
  ASSERT_TRUE(sol.converged);
  EXPECT_NEAR(sol.energy, -1.0 / 9.0, 1e-8);
}

TEST(RadialSolvers, ReportFailureWhenTheGridIsTooSmall) {
  // A 6s orbital of hydrogen does not fit in rmax = 5 bohr.
  const auto g = grid::RadialGrid::Logarithmic(-7.0, 0.01, 5.0, 1.0);
  const atom::PointNucleus nucleus(1.0);
  atom::Orbital o;
  o.n = 6;
  o.j = 0.5;
  const auto v = nucleus.Potential(g);
  EXPECT_FALSE(solvers::MakeSchrodingerSolver(g, nucleus)->Solve(o, v, -1.0 / 36).converged);
  EXPECT_FALSE(solvers::MakeDiracSolver(g, nucleus, constants::kSpeedOfLight)
                   ->Solve(o, v, -1.0 / 36).converged);
}

TEST(RadialGrid, IsOddAndRejectsBadInput) {
  const auto g = grid::RadialGrid::Logarithmic(-7.0, 0.0125, 60.0, 10.0);
  EXPECT_EQ(g.size() % 2, 1u);
  EXPECT_NEAR(g.r(0), std::exp(-7.0) / 10.0, 1e-18);
  EXPECT_THROW(grid::RadialGrid::Logarithmic(-7.0, 0.0, 60.0, 1.0), Error);
  const grid::RadialPoisson poisson(g);
  EXPECT_THROW(poisson.Solve(2, 2, g.Zeros()), Error);  // needs nst > k+1
}

}  // namespace
}  // namespace atomx
