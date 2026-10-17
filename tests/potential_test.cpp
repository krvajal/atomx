// Exact properties of self-consistent exact-exchange potentials.

#include <cmath>

#include <gtest/gtest.h>

#include "atomx/calc/calculation.hpp"
#include "atomx/core/constants.hpp"
#include "atomx/grid/poisson.hpp"
#include "atomx/xc/exchange_operator.hpp"
#include "tests/test_util.hpp"

namespace atomx {
namespace {

constexpr const char* kGrid = "xmin=-7, dx=0.0125, rmax=60, beta=0.3";

// For two electrons in one orbital, the exchange potential is exactly
// -1/2 of the Hartree potential (KLI and OEP are exact here).
void ExpectHalfHartree(const std::string& dft) {
  calc::Calculation calc(testing::InputFromNamelist(
      "zed=2., config='1s2', rel=0, dft='" + dft + "', " + kGrid));
  const auto res = calc.Run();
  ASSERT_TRUE(res.converged);
  const auto& g = calc.grid();
  const grid::RadialPoisson poisson(g);
  const Field vh = poisson.Solve(0, 2, res.density[0]);
  for (std::size_t i = 0; i < g.size(); i += 25) {
    if (g.r(i) > 10.0) break;
    EXPECT_NEAR(res.xc.potential[0][i], -0.5 * constants::kE2 * vh[i], 1e-7)
        << dft << " r=" << g.r(i);
  }
}

TEST(ExactExchangePotential, TwoElectronKliIsMinusHalfHartree) {
  ExpectHalfHartree("KLI");
}

TEST(ExactExchangePotential, TwoElectronOepIsMinusHalfHartree) {
  ExpectHalfHartree("OEP");
}

TEST(ExactExchangePotential, KliHomoConditionAndCoulombTail) {
  calc::Calculation calc(testing::InputFromNamelist(
      std::string("zed=10., config='default', rel=0, dft='KLI', ") + kGrid));
  const auto res = calc.Run();
  ASSERT_TRUE(res.converged);
  const auto& g = calc.grid();
  const auto& st = res.state;
  const Field& vx = res.xc.potential[0];

  // HOMO condition: <HOMO|Vx|HOMO> = <HOMO|K|HOMO>.
  std::size_t homo = 0;
  for (std::size_t k = 0; k < st.size(); ++k) {
    if (st.energy(k) > st.energy(homo)) homo = k;
  }
  const grid::RadialPoisson poisson(g);
  const auto op = xc::MakeNonRelativisticExchange(poisson);
  const Spinor k_homo = op->Apply(st, homo);
  const Field dens = st.OrbitalDensity(homo);
  Field w(g.size());
  for (std::size_t i = 0; i < w.size(); ++i) w[i] = dens[i] * vx[i];
  EXPECT_NEAR(g.Integrate(w, st.DensityPower(homo)),
              xc::ExchangeOperator::Expectation(st, homo, k_homo), 1e-8);

  // Asymptotics: Vx(r) -> -e^2 / r, approached from below with a 1/r^3
  // correction, so r Vx + e^2 shrinks monotonically and is ~1e-3 at 40 bohr.
  double previous = -1.0;
  for (std::size_t i = 0; i < g.size(); ++i) {
    if (g.r(i) < 8.0 || g.r(i) > 40.0) continue;
    const double gap = std::abs(vx[i] * g.r(i) + constants::kE2);
    if (previous >= 0.0) EXPECT_LE(gap, previous) << "r=" << g.r(i);
    previous = gap;
    EXPECT_LT(gap, 1.5 / (g.r(i) * g.r(i))) << "r=" << g.r(i);
  }
  EXPECT_LT(previous, 1e-3);
}

TEST(ExactExchangePotential, OepIsBelowKliAndSatisfiesVirial) {
  const std::string body = std::string("zed=12., config='default', rel=0, ") + kGrid;
  calc::Calculation kli(testing::InputFromNamelist(body + ", dft='KLI'"));
  calc::Calculation oep(testing::InputFromNamelist(body + ", dft='OEP'"));
  const auto rk = kli.Run();
  const auto ro = oep.Run();
  ASSERT_TRUE(rk.converged && ro.converged);
  // The OEP is the variational optimum among local potentials.
  EXPECT_LT(ro.energies.total, rk.energies.total);
  EXPECT_NEAR(0.5 * (rk.energies.total - ro.energies.total) * 1e3, 0.87, 0.01);
  EXPECT_NEAR(0.5 * ro.energies.VirialDiscrepancy() * 1e3, 0.0, 0.05);
}

}  // namespace
}  // namespace atomx
