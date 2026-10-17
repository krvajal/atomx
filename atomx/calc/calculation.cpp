#include "atomx/calc/calculation.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "atomx/core/constants.hpp"
#include "atomx/xc/kli.hpp"
#include "atomx/xc/lda.hpp"
#include "atomx/xc/oep.hpp"

namespace atomx::calc {

namespace {

std::unique_ptr<atom::Nucleus> BuildNucleus(const io::CalculationInput& in) {
  if (in.nucleus == "point") return std::make_unique<atom::PointNucleus>(in.z);
  if (in.nucleus == "sphere") {
    if (in.relativity != Relativity::kDirac) {
      throw Error("finite nucleus implemented only for rel=2");
    }
    const double radius = in.nuclear_radius > 0.0
                              ? in.nuclear_radius
                              : atom::UniformSphereNucleus::DefaultRadius(in.z);
    return std::make_unique<atom::UniformSphereNucleus>(in.z, radius);
  }
  throw Error("nucleus must be 'point' or 'sphere'");
}

atom::Configuration BuildConfiguration(const io::CalculationInput& in) {
  const atom::ConfigurationOptions opt{in.relativity, in.spin_polarized,
                                       in.distribution};
  return in.configuration == "default"
             ? atom::Configuration::Default(in.z, opt)
             : atom::Configuration::Parse(in.configuration, opt);
}

std::vector<std::string> SplitDft(const std::string& dft) {
  std::string s = dft;
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::toupper(c); });
  std::replace(s.begin(), s.end(), '+', ' ');
  std::replace(s.begin(), s.end(), '-', ' ');
  std::istringstream in(s);
  std::vector<std::string> parts;
  for (std::string t; in >> t;) parts.push_back(t);
  return parts;
}

}  // namespace

Calculation::Calculation(io::CalculationInput input)
    : input_(std::move(input)),
      grid_(grid::RadialGrid::Logarithmic(input_.xmin, input_.dx, input_.rmax,
                                          input_.z)),
      poisson_(grid_),
      nucleus_(BuildNucleus(input_)),
      configuration_(BuildConfiguration(input_)),
      speed_of_light_(input_.speed_of_light > 0.0 ? input_.speed_of_light
                                                  : constants::kSpeedOfLight),
      linear_solver_(linalg::MakeLinearSolver(input_.linear_solver)) {
  if (input_.relativity == Relativity::kDirac) {
    exchange_ = xc::MakeDiracExchange(poisson_);
    radial_solver_ = solvers::MakeDiracSolver(grid_, *nucleus_, speed_of_light_);
  } else {
    exchange_ = xc::MakeNonRelativisticExchange(poisson_);
    radial_solver_ = solvers::MakeSchrodingerSolver(grid_, *nucleus_);
  }
  xc_model_ = BuildXcModel(input_.dft);
}

std::unique_ptr<xc::XcModel> Calculation::BuildXcModel(const std::string& dft) {
  auto model = std::make_unique<xc::CompositeXcModel>();
  for (const auto& part : SplitDft(dft)) {
    if (part == "PZ" || part == "LDA") {
      // ld1.x's 'PZ' and 'LDA' both mean Slater exchange + PZ correlation;
      // after an orbital exchange ("KLI+PZ") only the correlation is added.
      if (part == "LDA" || SplitDft(dft).size() == 1) {
        model->Add(std::make_unique<xc::SlaterExchange>());
      }
      model->Add(std::make_unique<xc::PerdewZungerCorrelation>());
    } else if (part == "SLA") {
      model->Add(std::make_unique<xc::SlaterExchange>());
    } else if (part == "KLI") {
      model->Add(std::make_unique<xc::KliExchange>(*exchange_, *linear_solver_));
    } else if (part == "OEP") {
      if (input_.relativity != Relativity::kNone || !nucleus_->IsPointLike()) {
        throw Error("OEP is implemented only for rel=0 with a point nucleus");
      }
      response_ = std::make_unique<solvers::OrbitalResponse>(
          grid_, -nucleus_->charge() * constants::kE2);
      model->Add(std::make_unique<xc::OepExchange>(*exchange_, *response_,
                                                   *linear_solver_));
    } else {
      throw Error("unknown functional '" + part + "' in dft='" + dft + "'");
    }
  }
  return model;
}

scf::ScfResult Calculation::Run(scf::ScfDriver::Observer observer) {
  scf::ScfDriver driver(grid_, *nucleus_, *radial_solver_, *xc_model_,
                        scf::MakeAndersonMixer(input_.beta, input_.tr2),
                        scf::ScfOptions{input_.max_iterations});
  driver.SetObserver(std::move(observer));
  return driver.Run(configuration_, input_.relativity);
}

}  // namespace atomx::calc
