#include "atomx/report/report.hpp"

#include <cstdio>
#include <string>

#include "atomx/atom/periodic_table.hpp"
#include "atomx/core/constants.hpp"

namespace atomx::report {

namespace {

std::string Format(const char* fmt, double a, double b, double c) {
  char buf[160];
  std::snprintf(buf, sizeof buf, fmt, a, b, c);
  return buf;
}

void EnergyLine(std::ostream& out, const char* label, double ry) {
  char buf[160];
  std::snprintf(buf, sizeof buf, "     %-5s= %15.6f Ry, %15.6f Ha, %15.6f eV\n",
                label, ry, 0.5 * ry, ry * constants::kRydbergToEv);
  out << buf;
}

const char* RelativityName(Relativity r) {
  switch (r) {
    case Relativity::kNone: return "non-relativistic";
    case Relativity::kScalar: return "scalar-relativistic";
    case Relativity::kDirac: return "Dirac";
  }
  return "?";
}

}  // namespace

void WriteReport(std::ostream& out, const calc::Calculation& calc,
                 const scf::ScfResult& res) {
  const auto& in = calc.input();
  const auto& grid = calc.grid();
  char buf[200];

  out << "\n     atomx: all-electron atom\n\n";
  std::snprintf(buf, sizeof buf, "     atom %2s      Z = %3d   %s\n",
                atom::ElementSymbol(in.z).c_str(), in.z,
                RelativityName(in.relativity));
  out << buf;
  out << "     Exchange-correlation= " << calc.xc_model().Name() << "\n";
  out << "     nucleus: " << calc.nucleus().Describe() << "\n";
  std::snprintf(buf, sizeof buf,
                "     mesh =%5zu r(mesh) =%10.5f a.u. xmin =%6.2f dx =%8.5f\n",
                grid.size(), grid.r(grid.size() - 1), grid.xmin(), grid.dx());
  out << buf;
  if (in.relativity == Relativity::kDirac) {
    std::snprintf(buf, sizeof buf, "     c = %16.12f a.u.\n", calc.speed_of_light());
    out << buf;
  }

  const auto& st = res.state;
  const bool dirac = in.relativity == Relativity::kDirac;
  out << (dirac ? "\n     n l j   nl" : "\n     n l     nl")
      << "                  e(Ry)          e(Ha)          e(eV)\n";
  for (std::size_t k = 0; k < st.size(); ++k) {
    const auto& o = st.orbitals()[k];
    if (!o.IsActive()) continue;
    const double e = st.energy(k) - res.xc.eigenvalue_shift[o.spin];
    std::string head;
    if (dirac) {
      std::snprintf(buf, sizeof buf, "     %d %d %3.1f %2s %d(%6.2f)", o.n, o.l,
                    o.j, o.Label().c_str(), o.spin + 1, o.occupation);
    } else {
      std::snprintf(buf, sizeof buf, "     %d %d     %2s %d(%6.2f)", o.n, o.l,
                    o.Label().c_str(), o.spin + 1, o.occupation);
    }
    out << buf
        << Format("%15.4f%15.4f%15.4f\n", e, 0.5 * e,
                  e * constants::kRydbergToEv);
  }

  std::snprintf(buf, sizeof buf,
                "\n     final scf error: %8.1E reached in %3d iterations%s\n\n",
                res.residual, res.iterations,
                res.converged ? "" : "  (NOT CONVERGED)");
  out << buf;
  const auto& e = res.energies;
  EnergyLine(out, "Etot", e.total);
  out << "\n";
  EnergyLine(out, "Ekin", e.kinetic);
  EnergyLine(out, "Encl", e.nuclear);
  EnergyLine(out, "Eh", e.hartree);
  EnergyLine(out, "Exc", e.xc);
  std::snprintf(buf, sizeof buf,
                "\n     virial discrepancy -(Ekin+Etot) = %12.6f mHa\n",
                0.5 * e.VirialDiscrepancy() * 1000.0);
  out << buf;
}

}  // namespace atomx::report
