#ifndef BRANCH_BOUND_HPP
#define BRANCH_BOUND_HPP

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <llb.hpp>
#include <map>
#include <mokp.hpp>
#include <rwd.hpp>
#include <statistics.hpp>
#include <string>
#include <structs.hpp>
#include <subgraphs.hpp>
#include <time.hpp>
#include <tuple>
#include <utils.hpp>
#include <vector>

bool DEBUG = false;

std::vector<Statistics> HPS_BB(MOKP& instance, const int32_t K, const int32_t J, const double time_limit = 3600, const int32_t n_threads = 1) {
  if (DEBUG) {
    std::cout << "Number of objectives: " << instance.M << std::endl;
    std::cout << "Number of items: " << instance.N << std::endl;
    std::cout << "Number of solutions to return: " << J << std::endl;
  }

  // Set the number of threads
  omp_set_num_threads(n_threads);

  std::vector<Solution> solutions;
  std::vector<RefPoint> local_bounds = {RefPoint::origin(instance.M, instance.reference_point)};

  // Cache for the subproblems
  std::map<std::string, MOKPSolution> subproblems_cache;

  // Statistics
  std::vector<std::tuple<int32_t, int64_t, double>> statistics;
  const std::vector<int32_t> J_values = {5, 10, 15, 20, 50, 100};
  std::vector<Statistics> J_statistics;

  // Counters
  int64_t n_iterations = 0;
  int64_t n_subproblems = 0;
  int64_t n_cache_hits = 0;
  double subgraphs_time = 0.0;
  auto total_time = Timer(time_limit);

  while (!local_bounds.empty() && static_cast<int>(solutions.size()) < J) {
    // Check if the time limit has been reached
    if (total_time.finished()) {
      break;
    }
    // Save the statistics for the current size of the solutions if it is in J_values
    if (std::find(J_values.begin(), J_values.end(), static_cast<int>(solutions.size())) != J_values.end()) {
      J_statistics.push_back(Statistics(instance, solutions, statistics, total_time.elapsed(), n_iterations,
                                        n_subproblems, n_cache_hits, subgraphs_time, total_time.elapsed() - subgraphs_time));
    }

    if (DEBUG) {
      std::cout << "Iteration: " << n_iterations << std::endl;
      std::cout << "Current bounds:" << std::endl;
      for (const auto& local_bound : local_bounds) {
        std::cout << " - " << local_bound.to_string() << std::endl;
      }
    }

    MOKPSolution global_solution = MOKPSolution();
    for (int k = 1; k <= std::min(K, static_cast<int>(local_bounds.size())); k++) {
      // Check if the time limit has been reached
      if (total_time.finished()) {
        break;
      }

      if (DEBUG) {
        std::cout << "Subset size: " << k << std::endl;
      }

      auto subsets_time = Timer();
      std::vector<std::vector<int32_t>> subsets = Subgraphs::get_subgraphs(local_bounds, k, 120);
      subgraphs_time += subsets_time.elapsed();

      if (DEBUG) {
        std::cout << "Number of subsets: " << subsets.size() << std::endl;
      }

      MOKPSolution k_solution = MOKPSolution();

      std::mutex k_solution_mutex;  // Protect k_solution access

#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
      for (const auto& subset : subsets) {
        // Check if the time limit has been reached
        if (total_time.finished()) {
          continue; // Cannot use break in parallel for
        }

        if (DEBUG) {
          std::cout << "Subset: ";
          for (const auto& idx : subset) {
            std::cout << idx << " ";
          }
          std::cout << std::endl;
        }

        RefPoint union_point("up", std::vector<int64_t>(instance.M, 0));
        std::vector<RefPoint> ref_points;
        std::string subset_key = "-";
        for (const auto& ref_id : subset) {
          std::transform(local_bounds[ref_id].coordinates.begin(), local_bounds[ref_id].coordinates.end(),
                         union_point.coordinates.begin(), union_point.coordinates.begin(),
                         [](int64_t a, int64_t b) { return std::max(a, b); });
          ref_points.push_back(local_bounds[ref_id]);
          subset_key += local_bounds[ref_id].id + "-";
        }

        if (DEBUG) {
          std::cout << "Union point: " << union_point.to_string() << std::endl;
          std::cout << "Subset key: " << subset_key << std::endl;
          std::cout << "Ref points: " << ref_points.size() << std::endl;
          for (const auto& ref_point : ref_points) {
            std::cout << ref_point.to_string() << std::endl;
          }
        }

        MOKPSolution solution_found = MOKPSolution();
        // Protected cache access
        {
          std::lock_guard<std::mutex> lock(k_solution_mutex);
          if (subproblems_cache.find(subset_key) != subproblems_cache.end()) {
            solution_found = subproblems_cache[subset_key];
            n_cache_hits++;
          }
        }

        if (solution_found.hypervolume <= 0) {
          auto solution_time = Timer();
          const MOKPRecursiveState& state = instance.compute_solution(union_point, ref_points, total_time);
          solution_found = state.solution;
          solution_found.used_reference_points = ref_points;
          // Second critical section: update cache and stats
          std::lock_guard<std::mutex> lock(k_solution_mutex);
          n_subproblems++;
          subproblems_cache[subset_key] = solution_found;
          statistics.push_back(std::make_tuple(k, state.n_recursions, solution_time.elapsed()));
        }

        if (solution_found.hypervolume > 0) {
          std::lock_guard<std::mutex> lock(k_solution_mutex);
          if (DEBUG) {
            std::cout << "Found solution:" << solution_found.to_string() << std::endl;
          }
          if (solution_found.hypervolume > k_solution.hypervolume) {
            if (DEBUG) {
              std::cout << "Updating K solution..." << std::endl;
            }
            k_solution = solution_found;
          }
        } else {
          if (DEBUG) std::cout << "No solution found!" << std::endl;
        }
      }

      if (k_solution.hypervolume <= 0) {
        if (DEBUG) {
          std::cout << "Exiting k subsets..." << std::endl;
        }
        break;
      }
      if (k_solution.hypervolume > global_solution.hypervolume) {
        if (DEBUG) {
          std::cout << "Updating global solution..." << std::endl;
        }
        global_solution = k_solution;
      }
    }

    if (global_solution.hypervolume <= 0) {
      if (DEBUG) {
        std::cout << "Exiting main loop..." << std::endl;
      }
      break;
    }

    Solution solution = Solution("z" + std::to_string(solutions.size()), global_solution.hypervolume, 0, global_solution.values);
    solutions.push_back(solution);
    if (DEBUG) {
      std::cout << "Best solution found:" << solution.to_string() << std::endl;
    }

    // Update the local bounds
    LLB::update_bounds(local_bounds, solution);
    n_iterations++;

    if (DEBUG) {
      std::cout << "Local bounds: " << local_bounds.size() << std::endl;
      for (const auto& ref_point : local_bounds) {
        std::cout << ref_point.to_string() << std::endl;
      }
      std::cout << "Solutions: " << solutions.size() << std::endl;
      for (const auto& sol : solutions) {
        std::cout << sol.to_string() << std::endl;
      }
    }
  }

  double end_time = total_time.elapsed();
  double subproblems_time = end_time - subgraphs_time;

  // Save the statistics for the last size of the solutions (the J value received)
  J_statistics.push_back(Statistics(instance, solutions, statistics, end_time, n_iterations, n_subproblems, n_cache_hits, subgraphs_time, subproblems_time));

  if (DEBUG) {
    std::cout << "HPS Branch and Bound finished!" << std::endl;
    std::cout << "Total time: " << end_time << " seconds" << std::endl;
    std::cout << "Number of iterations: " << n_iterations << std::endl;
    std::cout << "Number of subproblems: " << n_subproblems << std::endl;
    std::cout << "Number of cache hits: " << n_cache_hits << std::endl;
    std::cout << "Subgraphs time: " << subgraphs_time << " seconds" << std::endl;
    std::cout << "Subproblems time: " << subproblems_time << " seconds" << std::endl;
    std::cout << "Total solutions: " << solutions.size() << std::endl;
  }

  return J_statistics;
}

#endif  // BRANCH_BOUND_HPP
