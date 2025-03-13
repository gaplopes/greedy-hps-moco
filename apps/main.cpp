#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <branch_bound.hpp>
#include <ilp.hpp>
#include <iostream>
#include <simulate.hpp>

// This function compares the solutions obtained using the branch-and-bound algorithm and the greedy algorithm
void compare_algorithms(MOKP& problem, const int32_t K, const int32_t J, const double time_limit) {
  std::vector<Statistics> stats_bb_J = HPS_BB(problem, K, J, time_limit);
  for (auto stats : stats_bb_J) {
    if (stats.solutions.size() == 0) {
      continue;
    }
    std::cout << "J=" << stats.solutions.size() << std::endl;
    std::cout << "BB" << std::endl;
    std::cout << stats.to_string() << std::endl;
    Statistics stats_sim_J = HPS_SIM(problem, stats.solutions.size());
    std::cout << "SIM_J" << std::endl;
    std::cout << stats_sim_J.to_string() << std::endl;
    Statistics compare(problem);
    std::cout << "BBvsSIM_J" << std::endl;
    compare.compare_solution_sets(stats.solutions, stats_sim_J.solutions);
  }
}

int main(int argc, char* argv[]) {
  CLI::App app{"Hypervolume Polychotomic Scheme for MOCO"};

  std::string instance_path;
  int32_t K, J;
  bool ilp = false, bb = false, greedy = false;
  double time_limit = 0.0;
  int32_t n_threads = 1;
  bool detailed_output = false;

  app.add_option("--input-file", instance_path, "Input file with the problem instance")
      ->required()
      ->check(CLI::ExistingFile);  // Ensure the file exists

  app.add_option("--K", K, "Number of reference points to use")
      ->default_val(3)
      ->check(CLI::PositiveNumber);  // Ensure criteria_limit is a positive number

  app.add_option("--J", J, "Number of solutions to obtain")
      ->default_val(-1)
      ->check(CLI::PositiveNumber);  // Ensure criteria_limit is a positive number

  app.add_option("--time-limit", time_limit, "Time limit in seconds")
      ->default_val(3600.0)
      ->check(CLI::PositiveNumber);  // Ensure criteria_limit is a positive number

  app.add_option("--n-threads", n_threads, "Number of threads to use")
      ->default_val(1)
      ->check(CLI::PositiveNumber);  // Ensure criteria_limit is a positive number

  app.add_option("--ilp", ilp, "Use the ILP algorithm")
      ->default_val(false);

  app.add_option("--bb", bb, "Use the branch-and-bound algorithm")
      ->default_val(true);

  app.add_option("--greedy", greedy, "Use the greedy algorithm")
      ->default_val(false);

  app.add_option("--detailed-output", detailed_output, "Print detailed output")
      ->default_val(false);

  CLI11_PARSE(app, argc, argv);

  std::cerr << "MOCO-HPS" << std::endl;
  std::cerr << "- Input file: " << instance_path << std::endl;
  std::cerr << "- K: " << K << std::endl;
  std::cerr << "- J: " << J << std::endl;
  std::cerr << "- Time limit: " << time_limit << std::endl;
  std::cerr << "- Number of threads: " << n_threads << std::endl;
  std::cerr << "- ILP: " << ilp << std::endl;
  std::cerr << "- BB: " << bb << std::endl;
  std::cerr << "- Greedy: " << greedy << std::endl;

  MOKP problem = MOKP::from_stream(std::ifstream(instance_path));

  if (bb) {
    // Solve the problem using the branch-and-bound algorithm
    std::vector<Statistics> stats_J = HPS_BB(problem, K, J, time_limit);
    Statistics stats = stats_J.back();
    std::cout << stats.to_string(detailed_output) << std::endl;
  } else if (greedy) {
    // Solve the problem using the greedy algorithm
    if (J == 0) {
      J = problem.nondominated_set.size();
    }
    Statistics stats = HPS_SIM(problem, J);
    std::cout << stats.to_string(detailed_output) << std::endl;
  } else {
    // Solve the problem using the ILP algorithm
    Statistics stats = HPS_ILP(problem, K, J);
    std::cout << stats.to_string(detailed_output) << std::endl;
  }

  return 0;
}

// Test the application with the following command:
// ./moco-hps --input-file={path} --J=10 --K=2
