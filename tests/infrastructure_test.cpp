// Mixer, input parsing and configuration rules.

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "atomx/atom/configuration.hpp"
#include "atomx/io/input.hpp"
#include "atomx/io/namelist.hpp"
#include "atomx/scf/mixer.hpp"
#include "tests/test_util.hpp"

namespace atomx {
namespace {

TEST(AndersonMixer, FindsFixedPointOfLinearMap) {
  // output = A x + b with a contracting, non-symmetric A.
  const std::vector<std::vector<double>> a = {
      {0.5, 0.2, 0.0}, {-0.1, 0.3, 0.1}, {0.0, 0.4, -0.6}};
  const std::vector<double> b = {1.0, -2.0, 0.5};
  auto mixer = scf::MakeAndersonMixer(0.3, 1e-24);
  std::vector<double> x(3, 0.0);
  bool converged = false;
  int iter = 1;
  for (; iter <= 60 && !converged; ++iter) {
    std::vector<double> out(3);
    for (int i = 0; i < 3; ++i) {
      out[i] = b[i];
      for (int j = 0; j < 3; ++j) out[i] += a[i][j] * x[j];
    }
    converged = mixer->Mix(out, x, iter).converged;
  }
  ASSERT_TRUE(converged);
  EXPECT_LT(iter, 30);  // much faster than linear mixing with beta=0.3
  for (int i = 0; i < 3; ++i) {
    double r = b[i] - x[i];
    for (int j = 0; j < 3; ++j) r += a[i][j] * x[j];
    EXPECT_NEAR(r, 0.0, 1e-11);
  }
}

TEST(Namelist, QuotesCommentsAndFortranReals) {
  const auto nl = io::Namelist::FromString(
      "junk before\n&INPUT\n"
      "  zed = 10., config='[He] 2s2, 2p6'  ! comment, with = sign\n"
      "  title=\"a/b\", beta=0.3d0\n"
      "/\n&other\n x=1\n/\n",
      "input");
  EXPECT_DOUBLE_EQ(nl.GetDouble("ZED", 0.0), 10.0);
  EXPECT_EQ(nl.GetString("config", ""), "[He] 2s2, 2p6");
  EXPECT_DOUBLE_EQ(nl.GetDouble("beta", 0.0), 0.3);
  EXPECT_EQ(nl.GetInt("missing", 7), 7);
  const auto unused = nl.Unused();
  ASSERT_EQ(unused.size(), 1u);
  EXPECT_EQ(unused.begin()->first, "title");
  EXPECT_EQ(unused.begin()->second, "a/b");
}

TEST(Namelist, Errors) {
  EXPECT_THROW(io::Namelist::FromString("&other\n/\n", "input"), Error);
  EXPECT_THROW(io::Namelist::FromString("&input\n zed=1\n", "input"), Error);
  const auto nl = io::Namelist::FromString("&input\n zed=abc\n/\n", "input");
  EXPECT_THROW(nl.GetDouble("zed", 0.0), Error);
}

TEST(CalculationInput, Ld1DefaultsAndChecks) {
  auto in = testing::InputFromNamelist("atom='Hg', rel=2, dft='KLI'");
  EXPECT_EQ(in.z, 80);
  EXPECT_EQ(in.relativity, Relativity::kDirac);
  EXPECT_DOUBLE_EQ(in.xmin, -8.0);  // ld1.x default for relativistic runs
  EXPECT_DOUBLE_EQ(in.dx, 0.008);
  EXPECT_EQ(in.configuration, "default");

  // rel defaults to scalar-relativistic above Z=18, which is not available
  EXPECT_THROW(testing::InputFromNamelist("zed=80."), Error);
  EXPECT_NO_THROW(testing::InputFromNamelist("zed=18."));
  EXPECT_THROW(testing::InputFromNamelist("zed=10., rel=1"), Error);
  EXPECT_THROW(testing::InputFromNamelist("zed=10.5"), Error);
  EXPECT_THROW(testing::InputFromNamelist("zed=10., atom='Ar'"), Error);
  EXPECT_THROW(testing::InputFromNamelist("zed=10., iswitch=3"), Error);
  EXPECT_THROW(testing::InputFromNamelist("zed=10., xmin=-1"), Error);
}

TEST(Configuration, CoresWithAndWithoutSpace) {
  const auto a = atom::Configuration::Parse("[Kr]5s2 4d10", {});
  const auto b = atom::Configuration::Parse("[Kr] 5s2 4d10", {});
  ASSERT_EQ(a.orbitals().size(), b.orbitals().size());
  EXPECT_DOUBLE_EQ(a.ElectronCount(), 48.0);
  EXPECT_DOUBLE_EQ(atom::Configuration::Parse("[Rn] 7s2", {}).ElectronCount(), 88.0);
}

TEST(Configuration, AverageSpinOrbitDistribution) {
  atom::ConfigurationOptions opt;
  opt.relativity = Relativity::kDirac;
  opt.distribution = atom::SpinOrbitDistribution::kAverage;
  const auto c = atom::Configuration::Parse("2p3", opt);
  ASSERT_EQ(c.orbitals().size(), 2u);
  EXPECT_DOUBLE_EQ(c.orbitals()[0].occupation, 1.0);  // 3 * 2/6
  EXPECT_DOUBLE_EQ(c.orbitals()[1].occupation, 2.0);  // 3 * 4/6
}

TEST(Configuration, ExplicitJAndSpinEntries) {
  atom::ConfigurationOptions dirac;
  dirac.relativity = Relativity::kDirac;
  const auto c = atom::Configuration::Parse("2p1 2p3", dirac);
  ASSERT_EQ(c.orbitals().size(), 2u);
  EXPECT_DOUBLE_EQ(c.orbitals()[0].j, 0.5);
  EXPECT_DOUBLE_EQ(c.orbitals()[1].j, 1.5);

  atom::ConfigurationOptions lsd;
  lsd.spin_polarized = true;
  const auto s = atom::Configuration::Parse("2p3 2p1", lsd);
  ASSERT_EQ(s.orbitals().size(), 2u);
  EXPECT_EQ(s.orbitals()[0].spin, 0);
  EXPECT_EQ(s.orbitals()[1].spin, 1);
}

TEST(Configuration, Errors) {
  EXPECT_THROW(atom::Configuration::Parse("", {}), Error);
  EXPECT_THROW(atom::Configuration::Parse("2p7", {}), Error);
  EXPECT_THROW(atom::Configuration::Parse("1p1", {}), Error);
  EXPECT_THROW(atom::Configuration::Parse("2x1", {}), Error);
  EXPECT_THROW(atom::Configuration::Parse("2p", {}), Error);
  EXPECT_THROW(atom::Configuration::Parse("[Zz] 1s1", {}), Error);
  EXPECT_THROW(atom::Configuration::Parse("1s1 1s1", {}), Error);
  atom::ConfigurationOptions dirac;
  dirac.relativity = Relativity::kDirac;
  EXPECT_THROW(atom::Configuration::Parse("2p3 2p3", dirac), Error);  // j=1/2 holds 2
  EXPECT_THROW(atom::Configuration::Parse("1s1 1s1", dirac), Error);
  dirac.spin_polarized = true;
  EXPECT_THROW(atom::Configuration::Parse("1s1", dirac), Error);
}

}  // namespace
}  // namespace atomx
