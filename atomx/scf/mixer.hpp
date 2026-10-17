#pragma once

#include <memory>
#include <vector>

namespace atomx::scf {

struct MixResult {
  bool converged = false;
  double residual = 0.0;  // mean squared difference output - input
};

// Builds the next SCF input from the current input and output potentials.
class Mixer {
 public:
  virtual ~Mixer() = default;
  // `input` is overwritten with the next input unless converged.
  // `iteration` starts at 1.
  virtual MixResult Mix(const std::vector<double>& output,
                        std::vector<double>& input, int iteration) = 0;
};

// Anderson mixing with up to two previous steps (ld1.x's dmixp, id=3).
std::unique_ptr<Mixer> MakeAndersonMixer(double beta, double threshold);

}  // namespace atomx::scf
