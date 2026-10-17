#include "atomx/report/data_writer.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <sstream>

#include "atomx/atom/periodic_table.hpp"
#include "atomx/core/types.hpp"

namespace atomx::report {

namespace {

std::ofstream Open(const std::string& path) {
  std::ofstream out(path);
  if (!out) throw Error("cannot write " + path);
  return out;
}

std::string Number(double x) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.12e", x);
  return buf;
}

std::string OrbitalName(const atom::Orbital& o, bool dirac, int nspin) {
  std::string name = o.Label();
  if (dirac) name += o.j < o.l ? "-" : "+";  // j = l -/+ 1/2
  if (nspin == 2) name += o.spin == 0 ? "up" : "dn";
  return name;
}

const char* RelativityName(Relativity r) {
  switch (r) {
    case Relativity::kNone: return "none";
    case Relativity::kScalar: return "scalar";
    case Relativity::kDirac: return "dirac";
  }
  return "?";
}

// Shared view of what is written.
struct Snapshot {
  const calc::Calculation& calc;
  const scf::ScfResult& res;
  bool dirac() const {
    return calc.input().relativity == Relativity::kDirac;
  }
  int nspin() const { return res.state.nspin(); }
  std::vector<std::size_t> active() const {
    std::vector<std::size_t> k;
    for (std::size_t i = 0; i < res.state.size(); ++i) {
      if (res.state.orbitals()[i].IsActive()) k.push_back(i);
    }
    return k;
  }
  double PrintedEnergy(std::size_t k) const {
    const auto& o = res.state.orbitals()[k];
    return res.state.energy(k) - res.xc.eigenvalue_shift[o.spin];
  }
  // Kohn-Sham potential that the final orbitals are eigenstates of.
  const Field& vks(int s) const { return res.state.potential()[s]; }
};

class ColumnWriter final : public ResultWriter {
 public:
  explicit ColumnWriter(std::string prefix) : prefix_(std::move(prefix)) {}

  std::vector<std::string> Write(const calc::Calculation& calc,
                                 const scf::ScfResult& res) const override {
    const Snapshot s{calc, res};
    const auto& g = calc.grid();
    const auto active = s.active();
    std::vector<std::string> files;

    auto header = [&](std::ofstream& out, const std::string& what) {
      out << "# atomx " << what << ": " << atom::ElementSymbol(calc.input().z)
          << " Z=" << calc.input().z << " dft=" << calc.xc_model().Name()
          << " rel=" << RelativityName(calc.input().relativity) << "\n";
    };

    {  // orbitals
      const std::string path = prefix_ + ".wfc";
      auto out = Open(path);
      header(out, "radial orbitals r*R(r), large component");
      out << "# energies (Ry):";
      for (auto k : active) out << " " << Number(s.PrintedEnergy(k));
      out << "\n#        r";
      for (auto k : active) {
        out << " " << OrbitalName(res.state.orbitals()[k], s.dirac(), s.nspin());
      }
      out << "\n";
      for (std::size_t i = 0; i < g.size(); ++i) {
        out << Number(g.r(i));
        for (auto k : active) out << " " << Number(res.state.spinor(k).large[i]);
        out << "\n";
      }
      files.push_back(path);
    }
    if (s.dirac()) {
      const std::string path = prefix_ + ".wfcs";
      auto out = Open(path);
      header(out, "radial orbitals, small component");
      out << "#        r";
      for (auto k : active) {
        out << " " << OrbitalName(res.state.orbitals()[k], true, s.nspin());
      }
      out << "\n";
      for (std::size_t i = 0; i < g.size(); ++i) {
        out << Number(g.r(i));
        for (auto k : active) out << " " << Number(res.state.spinor(k).small[i]);
        out << "\n";
      }
      files.push_back(path);
    }
    {  // densities
      const std::string path = prefix_ + ".rho";
      auto out = Open(path);
      header(out, "radial density 4 pi r^2 n(r)");
      out << "#        r" << (s.nspin() == 2 ? " up down" : " total") << "\n";
      for (std::size_t i = 0; i < g.size(); ++i) {
        out << Number(g.r(i));
        for (int sp = 0; sp < s.nspin(); ++sp) out << " " << Number(res.density[sp][i]);
        out << "\n";
      }
      files.push_back(path);
    }
    {  // potentials
      const std::string path = prefix_ + ".pot";
      auto out = Open(path);
      header(out, "potentials (Ry)");
      out << "#        r v_nucleus v_hartree";
      for (int sp = 0; sp < s.nspin(); ++sp) out << " v_xc" << sp + 1;
      for (int sp = 0; sp < s.nspin(); ++sp) out << " v_ks" << sp + 1;
      out << "\n";
      for (std::size_t i = 0; i < g.size(); ++i) {
        out << Number(g.r(i)) << " " << Number(res.v_nucleus[i]) << " "
            << Number(res.v_hartree[i]);
        for (int sp = 0; sp < s.nspin(); ++sp) out << " " << Number(res.xc.potential[sp][i]);
        for (int sp = 0; sp < s.nspin(); ++sp) out << " " << Number(s.vks(sp)[i]);
        out << "\n";
      }
      files.push_back(path);
    }
    {  // SCF history
      const std::string path = prefix_ + ".scf";
      auto out = Open(path);
      header(out, "SCF convergence");
      out << "# iteration residual\n";
      for (std::size_t i = 0; i < res.residual_history.size(); ++i) {
        out << i + 1 << " " << Number(res.residual_history[i]) << "\n";
      }
      files.push_back(path);
    }
    return files;
  }

 private:
  std::string prefix_;
};

// Minimal JSON emitter: objects, arrays of numbers and strings.
class Json {
 public:
  explicit Json(std::ostream& out) : out_(out) {}
  void Open(char c) { Separator(); out_ << c; first_ = true; }
  void Close(char c) { out_ << c; first_ = false; }
  void Key(const std::string& k) { Separator(); out_ << '"' << k << "\":"; first_ = true; }
  void Value(double x) { Separator(); out_ << Number(x); }
  void Value(int x) { Separator(); out_ << x; }
  void Value(bool x) { Separator(); out_ << (x ? "true" : "false"); }
  void Value(const std::string& x) { Separator(); out_ << '"' << x << '"'; }
  void Array(const Field& f) {
    Open('[');
    for (double x : f) Value(x);
    Close(']');
  }
  template <typename T>
  void Pair(const std::string& k, const T& v) { Key(k); Value(v); }

 private:
  void Separator() {
    if (!first_) out_ << ',';
    first_ = false;
  }
  std::ostream& out_;
  bool first_ = true;
};

class JsonWriter final : public ResultWriter {
 public:
  explicit JsonWriter(std::string prefix) : prefix_(std::move(prefix)) {}

  std::vector<std::string> Write(const calc::Calculation& calc,
                                 const scf::ScfResult& res) const override {
    const Snapshot s{calc, res};
    const auto& in = calc.input();
    const auto& g = calc.grid();
    const std::string path = prefix_ + ".json";
    auto out = Open(path);
    Json j(out);
    j.Open('{');
    j.Pair("program", std::string("atomx"));
    j.Pair("units", std::string("Rydberg atomic units; orbitals are r*R(r)"));
    j.Pair("element", atom::ElementSymbol(in.z));
    j.Pair("z", in.z);
    j.Pair("dft", calc.xc_model().Name());
    j.Pair("relativity", std::string(RelativityName(in.relativity)));
    j.Pair("nucleus", calc.nucleus().Describe());
    j.Pair("speed_of_light", calc.speed_of_light());
    j.Pair("nspin", s.nspin());
    j.Key("grid");
    j.Open('{');
    j.Pair("xmin", g.xmin());
    j.Pair("dx", g.dx());
    j.Key("r");
    j.Array(Field(g.r().begin(), g.r().end()));
    j.Close('}');

    j.Key("scf");
    j.Open('{');
    j.Pair("converged", res.converged);
    j.Pair("iterations", res.iterations);
    j.Key("residuals");
    j.Array(res.residual_history);
    j.Close('}');

    const auto& e = res.energies;
    j.Key("energies_ry");
    j.Open('{');
    j.Pair("total", e.total);
    j.Pair("kinetic", e.kinetic);
    j.Pair("nuclear", e.nuclear);
    j.Pair("hartree", e.hartree);
    j.Pair("xc", e.xc);
    j.Pair("virial_discrepancy", e.VirialDiscrepancy());
    j.Close('}');

    j.Key("orbitals");
    j.Open('[');
    for (auto k : s.active()) {
      const auto& o = res.state.orbitals()[k];
      j.Open('{');
      j.Pair("label", OrbitalName(o, s.dirac(), s.nspin()));
      j.Pair("n", o.n);
      j.Pair("l", o.l);
      j.Pair("j", o.j);
      j.Pair("spin", o.spin);
      j.Pair("occupation", o.occupation);
      j.Pair("energy", s.PrintedEnergy(k));
      j.Key("large");
      j.Array(res.state.spinor(k).large);
      if (s.dirac()) {
        j.Key("small");
        j.Array(res.state.spinor(k).small);
      }
      j.Close('}');
    }
    j.Close(']');

    j.Key("density");
    j.Open('[');
    for (int sp = 0; sp < s.nspin(); ++sp) j.Array(res.density[sp]);
    j.Close(']');

    j.Key("potentials_ry");
    j.Open('{');
    j.Key("nucleus");
    j.Array(res.v_nucleus);
    j.Key("hartree");
    j.Array(res.v_hartree);
    j.Key("xc");
    j.Open('[');
    for (int sp = 0; sp < s.nspin(); ++sp) j.Array(res.xc.potential[sp]);
    j.Close(']');
    j.Key("kohn_sham");
    j.Open('[');
    for (int sp = 0; sp < s.nspin(); ++sp) j.Array(s.vks(sp));
    j.Close(']');
    j.Close('}');
    j.Close('}');
    out << "\n";
    return {path};
  }

 private:
  std::string prefix_;
};

}  // namespace

std::unique_ptr<ResultWriter> MakeColumnWriter(std::string prefix) {
  return std::make_unique<ColumnWriter>(std::move(prefix));
}

std::unique_ptr<ResultWriter> MakeJsonWriter(std::string prefix) {
  return std::make_unique<JsonWriter>(std::move(prefix));
}

std::vector<std::unique_ptr<ResultWriter>> MakeWriters(const std::string& list,
                                                       const std::string& prefix) {
  std::string s = list;
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  std::replace(s.begin(), s.end(), '+', ' ');
  std::replace(s.begin(), s.end(), ',', ' ');
  std::istringstream in(s);
  std::vector<std::unique_ptr<ResultWriter>> writers;
  for (std::string t; in >> t;) {
    if (t == "dat") {
      writers.push_back(MakeColumnWriter(prefix));
    } else if (t == "json") {
      writers.push_back(MakeJsonWriter(prefix));
    } else if (t != "none") {
      throw Error("unknown output format '" + t + "' (use dat, json)");
    }
  }
  return writers;
}

}  // namespace atomx::report
