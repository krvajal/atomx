// Total energies against ld1.x (Quantum ESPRESSO 7.6), on the grid
// xmin=-7, dx=0.0125, rmax=60, beta=0.3. ld1.x prints 6 decimals in Ry, and
// atomx reproduces every one of them for these cases; 2e-6 Ry allows for
// rounding only.

#include <string>

#include <gtest/gtest.h>

#include "atomx/calc/calculation.hpp"
#include "tests/test_util.hpp"

namespace atomx {
namespace {

struct GoldenCase {
  const char* name;
  const char* input;  // &input body without the grid
  double etot_ry;     // ld1.x Etot
};

std::ostream& operator<<(std::ostream& os, const GoldenCase& c) {
  return os << c.name;
}

class GoldenLd1 : public ::testing::TestWithParam<GoldenCase> {};

TEST_P(GoldenLd1, TotalEnergyMatches) {
  const auto& c = GetParam();
  calc::Calculation calc(testing::InputFromNamelist(
      std::string(c.input) + ", xmin=-7, dx=0.0125, rmax=60, beta=0.3"));
  const auto res = calc.Run();
  ASSERT_TRUE(res.converged);
  EXPECT_NEAR(res.energies.total, c.etot_ry, 2e-6);
}

INSTANTIATE_TEST_SUITE_P(
    AllElectron, GoldenLd1,
    ::testing::Values(
        // LDA, LSDA, open shells
        GoldenCase{"ne_lda", "zed=10., config='default', rel=0, dft='PZ'", -256.454566},
        GoldenCase{"c_lda", "zed=6., config='default', rel=0, dft='PZ'", -74.848524},
        GoldenCase{"fe_lsda", "zed=26., config='default', rel=0, lsd=1, dft='PZ'", -2522.409254},
        // KLI, closed and open shells, spin-polarised
        GoldenCase{"c_kli", "zed=6., config='default', rel=0, dft='KLI'", -75.278163},
        GoldenCase{"o_kli", "zed=8., config='default', rel=0, dft='KLI'", -149.481263},
        GoldenCase{"ar_kli", "zed=18., config='default', rel=0, dft='KLI'", -1053.620959},
        GoldenCase{"au_kli", "zed=79., config='default', rel=0, dft='KLI'", -35730.598284},
        GoldenCase{"na_kli_lsd", "zed=11., config='default', rel=0, lsd=1, dft='KLI'", -323.711830},
        GoldenCase{"n_kli_lsd", "zed=7., config='default', rel=0, lsd=1, dft='KLI'", -108.806082},
        // OEP
        GoldenCase{"be_oep", "zed=4., config='default', rel=0, dft='OEP'", -29.144865},
        GoldenCase{"ne_oep", "zed=10., config='default', rel=0, dft='OEP'", -257.090831},
        GoldenCase{"n_oep_lsd", "zed=7., config='default', rel=0, lsd=1, dft='OEP'", -108.806794},
        // Dirac, point nucleus
        GoldenCase{"ar_dirac_lda", "zed=18., config='default', rel=2, dft='PZ'", -1055.618815},
        GoldenCase{"pt_dirac_lda", "zed=78., config='default', rel=2, dft='PZ'", -36872.917089},
        GoldenCase{"kr_rkli", "zed=36., config='default', rel=2, dft='KLI'", -5577.737750},
        GoldenCase{"u_rkli", "zed=92., config='default', rel=2, dft='KLI'", -56140.328851},
        // Dirac, finite nucleus
        GoldenCase{"ne_rkli_fn", "zed=10., config='default', rel=2, dft='KLI', nucleus='sphere'", -257.379151},
        GoldenCase{"kr_rkli_fn", "zed=36., config='default', rel=2, dft='KLI', nucleus='sphere'", -5577.689219},
        GoldenCase{"hg_lda_fn", "zed=80., config='default', rel=2, dft='PZ', nucleus='sphere'", -39293.171650}),
    [](const ::testing::TestParamInfo<GoldenCase>& info) {
      return std::string(info.param.name);
    });

}  // namespace
}  // namespace atomx
