#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ghps/hps/algorithm.hpp>
#include <ghps/hps/simulate.hpp>

using namespace moco::ghps;

int main(int argc, char *argv[]) {
  CLI::App app{"Greedy Hypervolume Polychotomic Scheme for MOCO"};

  std::string instance_path;
  int32_t K, J;
  double time_limit = 0.0;
  int32_t n_threads = 1;
  bool detailed_output = false;
  HPSMode HPS_mode = SIMULATION;

  app.add_option("--input-file", instance_path,
                 "Input file with the problem instance")
      ->required()
      ->check(CLI::ExistingFile); // Ensure the file exists

  app.add_option("--K", K, "Maximum reference points per scalarization")
      ->default_val(3)
      ->check(
          CLI::Range(1, 30)); // Reference subsets use signed integer masks

  app.add_option("--J", J, "Maximum positive-gain solutions to obtain")
      ->default_val(5)
      ->check(
          CLI::PositiveNumber);

  app.add_option("--time-limit", time_limit, "Time limit in seconds")
      ->default_val(3600.0)
      ->check(
          CLI::PositiveNumber);

  app.add_option("--n-threads", n_threads, "Number of OpenMP workers")
      ->default_val(1)
      ->check(
          CLI::PositiveNumber);

  app.add_option("--detailed-output", detailed_output, "Print detailed output")
      ->default_val(false);

  app.add_option(
         "--mode", HPS_mode,
         "HPS Mode: 0 = BB, 1 = QUADRATIC_HS, 2 = LOG_HS, 3 = HYBRID_HS, 4 = SIMULATION")
      ->default_val(SIMULATION)
      ->check(CLI::Range(0, 4));

  CLI11_PARSE(app, argc, argv);

  std::cerr << "Greedy-HPS-MOCO" << std::endl;
  std::cerr << "- Input file: " << instance_path << std::endl;
  std::cerr << "- K: " << K << std::endl;
  std::cerr << "- J: " << J << std::endl;
  std::cerr << "- Time limit: " << time_limit << std::endl;
  std::cerr << "- Number of threads: " << n_threads << std::endl;

  std::cerr << "- Detailed output: " << detailed_output << std::endl;
  std::cerr << "- HPS Mode: " << HPS_mode << std::endl;

  try {
    const MOKP problem = MOKP::from_stream(std::ifstream(instance_path));

    if (HPS_mode == SIMULATION) {
      std::cout << "[Simulation of the Greedy-HPS-MOCO]" << std::endl;
      // Solve the problem using the greedy algorithm
      if (problem.nondominated_set.empty()) throw std::invalid_argument("Simulation requires a supplied nondominated set");
      if (J > problem.n_nondominated_solutions) throw std::invalid_argument("Simulation J exceeds the supplied nondominated set size");
      Statistics stats = GHPS_SIM(problem, J);
      std::cout << stats.to_string(detailed_output) << std::endl;
      const auto filename = std::filesystem::path(instance_path).stem().string() + "_" + std::to_string(J) + ".seq";
      const auto dimension_folder = std::filesystem::path("output") / (std::to_string(problem.M) + "D");
      const auto folder = std::filesystem::is_directory(dimension_folder) ? dimension_folder : std::filesystem::path("output");
      const auto output_file = folder / filename;
      stats.save_sequence(output_file);
      std::cout << "Solutions saved to " << output_file.string() << std::endl;
    } else {
      std::cout << "[Greedy-HPS-MOCO with mode " << HPS_mode << "]" << std::endl;
      // Solve the problem using the selected scalarization backend
      std::vector<Statistics> stats_J = GHPS(problem, K, J, time_limit, n_threads, HPS_mode);
      for (const auto &stats : stats_J) {
        std::cout << "J=" << stats.solutions.size() << std::endl;
        std::cout << stats.to_string(detailed_output) << std::endl;
      }
      if (stats_J.back().status != SolverStatus::Optimal) return 1;
    }

    return 0;
  } catch (GRBException& e) {
    std::cerr << "Gurobi error: " << e.getMessage() << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 1;
}
