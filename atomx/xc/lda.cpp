#include "atomx/xc/lda.hpp"

#include <cmath>

#include "atomx/core/constants.hpp"

namespace atomx::xc {

namespace {

using constants::kE2;
using constants::kFourPi;

constexpr double kThreshold = 1.0e-10;      // XClib's rho_threshold_lda
constexpr double kPi34 = 0.6203504908994;   // (3/4pi)^(1/3)
constexpr double kThird = 1.0 / 3.0;
constexpr double kFourThirds = 4.0 / 3.0;

// Point values in Hartree: energy per electron and potential(s).
struct PointXc {
  double e = 0.0;
  double v_up = 0.0;
  double v_down = 0.0;
};

PointXc Slater(double n, double zeta, bool polarized) {
  constexpr double kAlpha = 2.0 / 3.0;
  PointXc p;
  if (!polarized) {
    constexpr double f = -0.687247939924714;  // -9/8 (3/2pi)^(2/3)
    const double rs = kPi34 / std::cbrt(n);
    p.e = f * kAlpha / rs;
    p.v_up = kFourThirds * f * kAlpha / rs;
    return p;
  }
  constexpr double f = -1.10783814957303361;
  const double up = f * kAlpha * std::cbrt((1.0 + zeta) * n);
  const double dw = f * kAlpha * std::cbrt((1.0 - zeta) * n);
  p.v_up = kFourThirds * up;
  p.v_down = kFourThirds * dw;
  p.e = 0.5 * ((1.0 + zeta) * up + (1.0 - zeta) * dw);
  return p;
}

struct PzParameters {
  double a, b, c, d, gc, b1, b2;
};
constexpr PzParameters kPzUnpolarized{0.0311, -0.048, 0.0020, -0.0116,
                                      -0.1423, 1.0529, 0.3334};
constexpr PzParameters kPzPolarized{0.01555, -0.0269, 0.0007, -0.0048,
                                    -0.0843, 1.3981, 0.2611};

void PerdewZunger(double rs, const PzParameters& p, double& ec, double& vc) {
  if (rs < 1.0) {
    const double lnrs = std::log(rs);
    ec = p.a * lnrs + p.b + p.c * rs * lnrs + p.d * rs;
    vc = p.a * lnrs + (p.b - p.a / 3.0) + 2.0 / 3.0 * p.c * rs * lnrs +
         (2.0 * p.d - p.c) / 3.0 * rs;
  } else {
    const double rs12 = std::sqrt(rs);
    const double ox = 1.0 + p.b1 * rs12 + p.b2 * rs;
    const double dox = 1.0 + 7.0 / 6.0 * p.b1 * rs12 + 4.0 / 3.0 * p.b2 * rs;
    ec = p.gc / ox;
    vc = ec * dox / ox;
  }
}

PointXc Pz(double n, double zeta, bool polarized) {
  PointXc out;
  const double rs = kPi34 / std::cbrt(n);
  double ecu, vcu;
  PerdewZunger(rs, kPzUnpolarized, ecu, vcu);
  if (!polarized) {
    out.e = ecu;
    out.v_up = vcu;
    return out;
  }
  double ecp, vcp;
  PerdewZunger(rs, kPzPolarized, ecp, vcp);
  const double denom = std::pow(2.0, kFourThirds) - 2.0;
  const double fz = (std::pow(1.0 + zeta, kFourThirds) +
                     std::pow(1.0 - zeta, kFourThirds) - 2.0) / denom;
  const double dfz = kFourThirds *
                     (std::pow(1.0 + zeta, kThird) - std::pow(1.0 - zeta, kThird)) /
                     denom;
  out.e = ecu + fz * (ecp - ecu);
  out.v_up = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (1.0 - zeta);
  out.v_down = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (-1.0 - zeta);
  return out;
}

// Evaluates a point functional on the radial grid; energy as \int e rho dr.
template <typename PointFunctional>
XcContribution EvaluateLocal(const atom::AtomState& state, const SpinField& rho,
                             PointFunctional functional) {
  const auto& grid = state.grid();
  const bool polarized = state.nspin() == 2;
  XcContribution out;
  out.potential = {grid.Zeros(), grid.Zeros()};
  Field energy_density = grid.Zeros();
  for (std::size_t i = 0; i < grid.size(); ++i) {
    const double norm = kFourPi * grid.r2()[i];
    const double up = rho[0][i] / norm;
    const double dw = polarized ? rho[1][i] / norm : 0.0;
    double n = polarized ? up + dw : std::abs(up);
    if (n <= kThreshold) continue;
    double zeta = 0.0;
    if (polarized) {
      zeta = (up - dw) / n;
      if (std::abs(zeta) > 1.0) zeta = zeta > 0 ? 1.0 : -1.0;
    }
    const PointXc p = functional(n, zeta, polarized);
    out.potential[0][i] = kE2 * p.v_up;
    if (polarized) out.potential[1][i] = kE2 * p.v_down;
    const double rho_total = rho[0][i] + (polarized ? rho[1][i] : 0.0);
    energy_density[i] = kE2 * p.e * rho_total;
  }
  out.energy = grid.Integrate(energy_density, 2);
  return out;
}

}  // namespace

XcContribution SlaterExchange::Evaluate(const atom::AtomState& state,
                                        const SpinField& rho) {
  return EvaluateLocal(state, rho, Slater);
}

XcContribution PerdewZungerCorrelation::Evaluate(const atom::AtomState& state,
                                                 const SpinField& rho) {
  return EvaluateLocal(state, rho, Pz);
}

}  // namespace atomx::xc
