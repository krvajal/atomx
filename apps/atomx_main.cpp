// atomx: reads an ld1.x-style &input namelist from stdin (or a file given as
// the first argument), runs the all-electron SCF and prints the results.

#include <fstream>
#include <iostream>

#include "atomx/calc/calculation.hpp"
#include "atomx/io/input.hpp"
#include "atomx/io/namelist.hpp"
#include "atomx/report/data_writer.hpp"
#include "atomx/report/report.hpp"

int main(int argc, char** argv) {
  try {
    atomx::io::Namelist nl;
    if (argc > 1) {
      std::ifstream file(argv[1]);
      if (!file) throw atomx::Error(std::string("cannot open ") + argv[1]);
      nl = atomx::io::Namelist::Read(file, "input");
    } else {
      nl = atomx::io::Namelist::Read(std::cin, "input");
    }
    auto input = atomx::io::CalculationInput::FromNamelist(nl);
    for (const auto& [key, value] : nl.Unused()) {
      std::cerr << "warning: ignoring unknown input '" << key << "'\n";
    }
    // Validate the output list before the (possibly long) run.
    const auto writers = atomx::report::MakeWriters(input.output, input.prefix);
    atomx::calc::Calculation calculation(std::move(input));
    const auto result = calculation.Run();
    atomx::report::WriteReport(std::cout, calculation, result);
    for (const auto& writer : writers) {
      for (const auto& path : writer->Write(calculation, result)) {
        std::cout << "     wrote " << path << "\n";
      }
    }
    return result.converged ? 0 : 2;
  } catch (const std::exception& e) {
    std::cerr << "atomx: error: " << e.what() << "\n";
    return 1;
  }
}
