#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <fstream>
#include <iostream>
#include "experiments.hpp"

using namespace moco::ghps;

int main(int argc, char *argv[]) {
  CLI::App app{"Greedy Hypervolume Polychotomic Scheme for the MOKP - Tests"};

  std::string instance_path;
  TestConfig config;

  app.add_option("--input-file", instance_path,
                 "Input file with the problem instance")
      ->required()
      ->check(CLI::ExistingFile);

  app.add_option("--seed", config.seed, "Random seed for the algorithms")
      ->default_val(42)
      ->check(CLI::PositiveNumber);

  app.add_flag("--milp", config.solvers.at(Solver::Milp), "Solve using MILP");
  app.add_flag("--miqp", config.solvers.at(Solver::Miqp), "Solve using MIQP");
  app.add_flag("--log", config.solvers.at(Solver::Log), "Solve using Log-Exp");
  app.add_flag("--greedy", config.solvers.at(Solver::Greedy), "Solve using Greedy");
  app.add_flag("--bb", config.solvers.at(Solver::BranchBound), "Solve using Branch and Bound");
      
  app.add_option("--threads", config.solver_threads, "Gurobi threads per solve")
      ->default_val(4)
      ->check(CLI::PositiveNumber);
      
  app.add_option("--refs", config.num_ref_points, "Number of reference points")
      ->default_val(30)
      ->check(CLI::PositiveNumber);

  app.add_option("--max-k", config.max_k, "Maximum number of reference points in combinations")
      ->default_val(5)
      ->check(CLI::PositiveNumber);

  app.add_option("--max-permutations", config.max_permutations, "Maximum permutations per k")
      ->default_val(30)
      ->check(CLI::PositiveNumber);

  CLI11_PARSE(app, argc, argv);

  std::cerr << "Running Tests" << std::endl;
  std::cerr << "- Input file: " << instance_path << std::endl;
  std::cerr << "- Seed: " << config.seed << std::endl;

  try {
    const MOKP problem = MOKP::from_stream(std::ifstream(instance_path));
    experiments(problem, config);

    return 0;
  } catch (GRBException& e) {
    std::cerr << "Gurobi error: " << e.getMessage() << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 1;
}
