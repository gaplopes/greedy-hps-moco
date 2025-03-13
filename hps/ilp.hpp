#ifndef ILP_HPP
#define ILP_HPP

#include <mokp.hpp>
#include <llb.hpp>
#include <structs.hpp>
#include <utils.hpp>
#include <vector>
#include <statistics.hpp>
// #include <gurobi.hpp>

Statistics HPS_ILP(const MOKP &instance, const int32_t &K = 3, const int32_t &MAX_TIME = 3600) {
  // Get the instance parameters
  int32_t DIM = instance.M;

  // Vector of local bounds
  std::vector<RefPoint> local_bounds = {RefPoint::origin(DIM, instance.reference_point)};

  // Vector of solutions
  std::vector<Solution> solutions;

  auto total_time = Timer();
  // HPS algorithm using ILP
  while ((int)local_bounds.size() > 0) {
    int64_t best_hypervolume = 0;
    Solution best_solution = Solution();
    for (int k = 1; k <= std::min(K, static_cast<int>(local_bounds.size())); k++) {
      int64_t best_k_hypervolume = 0;
      Solution best_k_solution;
      std::vector<std::vector<int32_t>> subsets = Utils::get_k_subsets(local_bounds.size(), k);
      for (const auto &subset : subsets) {
        std::vector<RefPoint> subset_bounds;
        for (auto idx : subset) subset_bounds.push_back(local_bounds[idx]);
        // Solve the hypervolume scalarization problem using the subset of k bounds
        bool finished = true;
        // auto gurobi_solution = solveILP(instance, subset_bounds, DIM, k, MAX_TIME, finished); // TODO: Install Gurobi and uncomment this line
        std::tuple<int64_t, std::string, std::vector<int64_t>> gurobi_solution;
        if (finished) {
          int64_t sol_hypervolume = std::get<0>(gurobi_solution);
          std::vector<int64_t> sol_coordinates = std::get<2>(gurobi_solution);
          if (sol_hypervolume > best_k_hypervolume) {
            best_k_hypervolume = sol_hypervolume;
            std::string best_k_solution_id = "z" + std::to_string(solutions.size());
            best_k_solution = Solution(best_k_solution_id, best_k_hypervolume, 0, sol_coordinates, subset_bounds);
          }
        }
      }
      if (best_k_hypervolume == 0) break;
      if (best_k_hypervolume > best_hypervolume) {
        best_hypervolume = best_k_hypervolume;
        best_solution = best_k_solution;
      }
    }
    if (best_hypervolume == 0) break;
    solutions.push_back(best_solution);
    // Update the local bounds
    LLB::update_bounds(local_bounds, best_solution);
  }
  auto end_time = total_time.elapsed();

  // Return the solutions found
  return Statistics(instance, solutions, {}, end_time, 0, 0, 0, 0, 0);
}

#endif  // ILP_HPP
