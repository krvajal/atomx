// The data files written for plotting: content and layout.

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "atomx/calc/calculation.hpp"
#include "atomx/report/data_writer.hpp"
#include "tests/test_util.hpp"

namespace atomx {
namespace {

std::string TempPrefix(const std::string& name) {
  const char* dir = std::getenv("TEST_TMPDIR");
  return std::string(dir ? dir : "/tmp") + "/" + name;
}

std::string Slurp(const std::string& path) {
  std::ifstream in(path);
  EXPECT_TRUE(in.good()) << path;
  std::stringstream s;
  s << in.rdbuf();
  return s.str();
}

TEST(DataWriter, ColumnFilesHaveOneRowPerGridPoint) {
  calc::Calculation calc(testing::InputFromNamelist(
      "zed=7., config='default', rel=0, lsd=1, dft='KLI', "
      "xmin=-7, dx=0.0125, rmax=60, beta=0.3"));
  const auto res = calc.Run();
  const std::string prefix = TempPrefix("n_kli");
  const auto files = report::MakeColumnWriter(prefix)->Write(calc, res);
  EXPECT_EQ(files.size(), 4u);  // no small components without Dirac

  // .wfc: 2 header lines + column names, then r and one column per active
  // orbital (N, spin-polarised: 1s, 2s, 2p up and down).
  std::istringstream wfc(Slurp(prefix + ".wfc"));
  std::string line;
  std::size_t rows = 0, cols = 0;
  while (std::getline(wfc, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream fields(line);
    std::size_t n = 0;
    for (std::string f; fields >> f;) ++n;
    if (rows == 0) cols = n;
    EXPECT_EQ(n, cols);
    ++rows;
  }
  EXPECT_EQ(rows, calc.grid().size());
  EXPECT_EQ(cols, 1 + res.state.size());

  // .pot: r, v_nucleus, v_hartree, v_xc up/down, v_ks up/down
  EXPECT_NE(Slurp(prefix + ".pot").find("v_xc2 v_ks1 v_ks2"), std::string::npos);
  // .scf: one line per iteration
  std::istringstream scf(Slurp(prefix + ".scf"));
  std::size_t iters = 0;
  while (std::getline(scf, line)) iters += (!line.empty() && line[0] != '#');
  EXPECT_EQ(iters, static_cast<std::size_t>(res.iterations));
}

TEST(DataWriter, JsonHasDiracSmallComponentsAndEnergies) {
  calc::Calculation calc(testing::InputFromNamelist(
      "zed=10., config='default', rel=2, dft='KLI', nucleus='sphere', "
      "xmin=-7, dx=0.0125, rmax=60, beta=0.3"));
  const auto res = calc.Run();
  const std::string prefix = TempPrefix("ne_rkli");
  report::MakeJsonWriter(prefix)->Write(calc, res);
  const std::string json = Slurp(prefix + ".json");
  EXPECT_EQ(json.front(), '{');
  EXPECT_EQ(json[json.size() - 2], '}');
  EXPECT_NE(json.find("\"relativity\":\"dirac\""), std::string::npos);
  EXPECT_NE(json.find("\"label\":\"2P-\""), std::string::npos);  // 2p1/2
  EXPECT_NE(json.find("\"small\":["), std::string::npos);
  EXPECT_NE(json.find("\"virial_discrepancy\""), std::string::npos);
  EXPECT_NE(json.find("uniformly charged sphere"), std::string::npos);
}

TEST(DataWriter, FormatListIsValidated) {
  EXPECT_EQ(report::MakeWriters("dat+json", "x").size(), 2u);
  EXPECT_EQ(report::MakeWriters("none", "x").size(), 0u);
  EXPECT_THROW(report::MakeWriters("hdf5", "x"), Error);
}

}  // namespace
}  // namespace atomx
