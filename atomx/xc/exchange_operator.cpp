#include "atomx/xc/exchange_operator.hpp"

#include <cmath>
#include <cstdlib>

#include "atomx/core/constants.hpp"
#include "atomx/math/wigner.hpp"

namespace atomx::xc {

using constants::kE2;

std::vector<Spinor> ExchangeOperator::ApplyAll(
    const atom::AtomState& state) const {
  std::vector<Spinor> out;
  out.reserve(state.size());
  for (std::size_t k = 0; k < state.size(); ++k) out.push_back(Apply(state, k));
  return out;
}

double ExchangeOperator::Expectation(const atom::AtomState& state,
                                     std::size_t k, const Spinor& k_psi) {
  const auto& psi = state.spinor(k);
  Field w(state.grid().size());
  for (std::size_t i = 0; i < w.size(); ++i) {
    w[i] = psi.large[i] * k_psi.large[i] + psi.small[i] * k_psi.small[i];
  }
  return state.grid().Integrate(w, state.DensityPower(k));
}

double ExchangeOperator::Energy(const atom::AtomState& state,
                                const std::vector<Spinor>& k_psi) {
  double e = 0.0;
  for (std::size_t k = 0; k < state.size(); ++k) {
    const double occ = state.orbitals()[k].occupation;
    if (occ > 0.0) e += 0.5 * occ * Expectation(state, k, k_psi[k]);
  }
  return e;
}

namespace {

// Adds fac * V_lam[psi_mu psi_nu] * psi_mu to `out`.
void AddTerm(const grid::RadialPoisson& poisson, int lam, int nst,
             const Field& pair_density, double fac, const Spinor& psi_mu,
             Spinor& out) {
  const Field v = poisson.Solve(lam, nst, pair_density);
  for (std::size_t i = 0; i < v.size(); ++i) {
    out.large[i] += fac * v[i] * psi_mu.large[i];
    out.small[i] += fac * v[i] * psi_mu.small[i];
  }
}

class NonRelativisticExchange final : public ExchangeOperator {
 public:
  explicit NonRelativisticExchange(const grid::RadialPoisson& poisson)
      : poisson_(poisson) {}

  Spinor Apply(const atom::AtomState& state, std::size_t nu) const override {
    const auto& orbitals = state.orbitals();
    const std::size_t mesh = state.grid().size();
    Spinor out{Field(mesh, 0.0), Field(mesh, 0.0)};
    const int l1 = orbitals[nu].l;
    const double half = 2.0 * l1 + 1.0;
    const auto& psi_nu = state.spinor(nu);
    for (std::size_t mu = 0; mu < state.size(); ++mu) {
      if (orbitals[mu].spin != orbitals[nu].spin) continue;
      // occupation per spin channel
      double ocs = orbitals[mu].occupation * 0.5 * state.nspin();
      double doc = 0.0;
      if (mu == nu && l1 != 0 && ocs > 0.0) {
        // spherical average of a partially filled shell
        const int i = static_cast<int>(ocs);
        doc = (i * (2.0 * ocs - i - 1.0) / (half - 1.0) - ocs * ocs / half) *
              half / ocs;
      }
      ocs += doc;
      if (ocs == 0.0 && doc == 0.0) continue;
      const int l2 = orbitals[mu].l;
      const auto& psi_mu = state.spinor(mu);
      Field pair(mesh);
      for (std::size_t i = 0; i < mesh; ++i) {
        pair[i] = psi_mu.large[i] * psi_nu.large[i];
      }
      for (int l3 = std::abs(l1 - l2); l3 <= l1 + l2; ++l3) {
        const double sss = math::LegendreTripleIntegral(l1, l2, l3);
        if (std::abs(sss) <= 1e-10) continue;
        AddTerm(poisson_, l3, l1 + l2 + 2, pair, -kE2 * ocs * sss / 2.0,
                psi_mu, out);
      }
      if (mu == nu && doc != 0.0) {
        // remove the spurious Hartree part of the average
        AddTerm(poisson_, 0, 2, pair, doc * kE2, psi_mu, out);
      }
    }
    return out;
  }

 private:
  const grid::RadialPoisson& poisson_;
};

class DiracExchange final : public ExchangeOperator {
 public:
  explicit DiracExchange(const grid::RadialPoisson& poisson)
      : poisson_(poisson) {}

  Spinor Apply(const atom::AtomState& state, std::size_t nu) const override {
    const auto& orbitals = state.orbitals();
    const std::size_t mesh = state.grid().size();
    Spinor out{Field(mesh, 0.0), Field(mesh, 0.0)};
    const auto& o_nu = orbitals[nu];
    const int two_j1 = static_cast<int>(std::lround(2.0 * o_nu.j));
    const double twoj1 = two_j1;
    const auto& psi_nu = state.spinor(nu);
    for (std::size_t mu = 0; mu < state.size(); ++mu) {
      const auto& o_mu = orbitals[mu];
      if (o_mu.occupation <= 0.0 && mu != nu) continue;
      const int two_j2 = static_cast<int>(std::lround(2.0 * o_mu.j));
      const auto& psi_mu = state.spinor(mu);
      Field pair(mesh);
      for (std::size_t i = 0; i < mesh; ++i) {
        pair[i] = psi_mu.large[i] * psi_nu.large[i] +
                  psi_mu.small[i] * psi_nu.small[i];
      }
      // G and F go roughly as r^(j+1/2) at small r
      const int nst = (two_j1 + two_j2) / 2 + 1;
      for (int lam = std::abs(two_j1 - two_j2) / 2; lam <= (two_j1 + two_j2) / 2;
           ++lam) {
        if ((o_nu.l + o_mu.l + lam) % 2 != 0) continue;
        const double cg =
            math::ClebschGordan(two_j1, 1, two_j2, -1, 2 * lam, 0);
        const double bb = cg * cg / (2.0 * lam + 1.0);
        double fac;
        if (mu != nu) {
          fac = -kE2 * o_mu.occupation * bb;
        } else {
          const double occ = o_nu.occupation;
          fac = (occ - 1.0) * (twoj1 + 1.0) / twoj1 * bb;
          if (lam == 0) fac += (twoj1 + 1.0 - occ) / twoj1;
          fac *= -kE2;
        }
        if (std::abs(fac) < 1e-12) continue;
        // the Poisson series needs nst > lam + 1
        AddTerm(poisson_, lam, std::max(nst, lam + 2), pair, fac, psi_mu, out);
      }
    }
    return out;
  }

 private:
  const grid::RadialPoisson& poisson_;
};

}  // namespace

std::unique_ptr<ExchangeOperator> MakeNonRelativisticExchange(
    const grid::RadialPoisson& poisson) {
  return std::make_unique<NonRelativisticExchange>(poisson);
}

std::unique_ptr<ExchangeOperator> MakeDiracExchange(
    const grid::RadialPoisson& poisson) {
  return std::make_unique<DiracExchange>(poisson);
}

}  // namespace atomx::xc
