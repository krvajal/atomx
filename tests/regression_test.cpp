// End-to-end SCF runs against published exact-exchange results.
//
// KLI and OEP: Ryabinkin, Kananenka, Staroverov, PRL 111, 013001 (2013).
// Relativistic KLI: Kreibich, Gross, Engel, PRA 57, 138 (1998), Table I.
// LDA: ld1.x (Quantum ESPRESSO 7.6) on the same grid.

#include <memory>

#include <gtest/gtest.h>

#include "atomx/atom/nucleus.hpp"
#include "atomx/calc/calculation.hpp"
#include "atomx/core/constants.hpp"
#include "atomx/grid/poisson.hpp"
#include "atomx/linalg/linear_solver.hpp"
#include "atomx/scf/scf_driver.hpp"
#include "atomx/solvers/radial_solver.hpp"
#include "atomx/xc/exchange_operator.hpp"
#include "atomx/xc/kli.hpp"

namespace atomx {
namespace {

io::CalculationInput Input(int z, const std::string& dft, int rel = 0,
                           bool lsd = false) {
  io::CalculationInput in;
  in.z = z;
  in.dft = dft;
  in.relativity = rel == 2 ? Relativity::kDirac : Relativity::kNone;
  in.spin_polarized = lsd;
  in.xmin = -7.0;
  in.dx = 0.0125;
  in.rmax = 60.0;
  in.beta = 0.3;
  return in;
}

double TotalHartree(io::CalculationInput in) {
  calc::Calculation calc(std::move(in));
  const auto r = calc.Run();
  EXPECT_TRUE(r.converged);
  return 0.5 * r.energies.total;
}

TEST(Regression, LdaNeonMatchesLd1) {
  EXPECT_NEAR(TotalHartree(Input(10, "PZ")), -128.227283, 1e-6);
}

TEST(Regression, KliNeonEnergyAndVirial) {
  calc::Calculation calc(Input(10, "KLI"));
  const auto r = calc.Run();
  ASSERT_TRUE(r.converged);
  EXPECT_NEAR(0.5 * r.energies.total, -128.54483, 1e-5);
  // KLI violates the Levy-Perdew virial relation: 155.62 mHa for Ne.
  EXPECT_NEAR(0.5 * r.energies.VirialDiscrepancy() * 1e3, 155.62, 0.01);
}

TEST(Regression, KliSpinPolarisedNitrogen) {
  EXPECT_NEAR(TotalHartree(Input(7, "KLI", 0, true)), -54.40304, 1e-5);
}

TEST(Regression, OepNeonSatisfiesVirial) {
  calc::Calculation calc(Input(10, "OEP"));
  const auto r = calc.Run();
  ASSERT_TRUE(r.converged);
  EXPECT_NEAR(0.5 * r.energies.total, -128.54541, 1e-5);
  EXPECT_NEAR(0.5 * r.energies.VirialDiscrepancy() * 1e3, 0.0, 0.1);
}

TEST(Regression, RelativisticKliFiniteNucleusKrypton) {
  auto in = Input(36, "KLI", 2);
  in.nucleus = "sphere";
  in.speed_of_light = 137.0359895;
  EXPECT_NEAR(TotalHartree(in), -2788.845, 1e-3);
}

TEST(Regression, RelativisticKliFiniteNucleusMercury) {
  auto in = Input(80, "KLI", 2);
  in.nucleus = "sphere";
  in.speed_of_light = 137.0359895;
  EXPECT_NEAR(TotalHartree(in), -19648.815, 1e-3);
}

// Components are plain interfaces: here KLI runs with a test double for the
// linear solver, wired by hand instead of through calc::Calculation.
class CountingSolver final : public linalg::LinearSolver {
 public:
  std::string name() const override { return "counting"; }
  std::vector<double> Solve(const linalg::DenseMatrix& a,
                            const std::vector<double>& b,
                            linalg::MatrixStructure s) const override {
    ++calls;
    return inner->Solve(a, b, s);
  }
  std::unique_ptr<linalg::LinearSolver> inner = linalg::MakeLapackSolver();
  mutable int calls = 0;
};

TEST(Composition, KliWithInjectedLinearSolver) {
  const auto grid = grid::RadialGrid::Logarithmic(-7.0, 0.0125, 60.0, 12.0);
  const grid::RadialPoisson poisson(grid);
  const atom::PointNucleus nucleus(12.0);
  const auto radial = solvers::MakeSchrodingerSolver(grid, nucleus);
  const auto exchange = xc::MakeNonRelativisticExchange(poisson);
  CountingSolver solver;
  xc::KliExchange kli(*exchange, solver);
  scf::ScfDriver driver(grid, nucleus, *radial, kli,
                        scf::MakeAndersonMixer(0.3, 1e-14));
  const auto r = driver.Run(atom::Configuration::Default(12, {}),
                            Relativity::kNone);
  ASSERT_TRUE(r.converged);
  EXPECT_NEAR(0.5 * r.energies.total, -199.61071, 1e-5);  // Mg, KLI
  EXPECT_EQ(solver.calls, r.iterations);  // one KLI solve per SCF step
}

}  // namespace
}  // namespace atomx
