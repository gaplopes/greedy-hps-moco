#ifndef SIMULATE_HPP
#define SIMULATE_HPP

#include <hvi.hpp>
#include <iostream>
#include <limits>
#include <llb.hpp>
#include <mokp.hpp>
#include <statistics.hpp>
#include <structs.hpp>
#include <utils.hpp>
#include <vector>

Statistics HPS_SIM(const MOKP& instance, int32_t J = -1) {
  // Verify if there are nondominated solutions
  if (instance.nondominated_set.empty()) {
    std::cout << "No nondominated solutions found." << std::endl;
    return Statistics(instance, {}, {}, 0, 0, 0, 0, 0, 0);
  }

  // Get the instance parameters
  int32_t M = instance.M;

  // Set J to the number of nondominated solutions if it is -1
  if (J == -1) {
    J = instance.nondominated_set.size();
  }

  std::vector<Solution> n_set = instance.nondominated_set;

  // Vector of local bounds
  RefPoint ref_point = RefPoint::origin(M, instance.reference_point);
  std::vector<RefPoint> local_bounds = {ref_point};

  // Vector of solutions
  std::vector<Solution> solutions;
  HypervolumeIndicator<int64_t, std::vector<int64_t>> hvc_space(ref_point.coordinates, true);

  auto total_time = Timer();
  // Simulate the HPS algorithm
  while (!n_set.empty() && static_cast<int>(solutions.size()) < J) {
    int64_t best_hypervolume = std::numeric_limits<int64_t>::min();
    Solution best_solution;  // Default constructor
    int32_t best_solution_index = -1;

    // For each solution, compute the hypervolume contribution using all dominated bounds
    for (size_t i = 0; i < n_set.size(); ++i) {
      Solution& sol = n_set[i];
      std::vector<RefPoint> dominated_bounds = LLB::get_dominated_bounds(local_bounds, sol);
      if (dominated_bounds.empty()) continue;

      int64_t sol_hypervolume = hvc_space.contribution(sol.coordinates);
      if (sol_hypervolume > best_hypervolume) {
        best_hypervolume = sol_hypervolume;
        best_solution = sol;
        best_solution.hypervolume = sol_hypervolume;
        best_solution_index = static_cast<int32_t>(i);
      }
    }

    if (best_solution_index == -1) break;  // There are no solutions to add

    // Store the results
    std::string id = "z" + std::to_string(solutions.size());
    std::vector<RefPoint> dominated_bounds = LLB::get_dominated_bounds(local_bounds, best_solution);
    Solution sol = Solution(id, best_hypervolume, 0, best_solution.coordinates, dominated_bounds);
    solutions.push_back(sol);
    hvc_space.insert(sol.coordinates);
    // Remove the best solution from the solutions vector
    n_set.erase(n_set.begin() + best_solution_index);
    // Update the local bounds
    LLB::update_bounds(local_bounds, best_solution);
  }
  auto end_time = total_time.elapsed();

  return Statistics(instance, solutions, {}, end_time, 0, 0, 0, 0, 0);
}

#endif  // SIMULATE_HPP
