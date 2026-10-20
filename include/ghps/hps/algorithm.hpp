#ifndef GHPS_HPS_ALGORITHM_HPP
#define GHPS_HPS_ALGORITHM_HPP

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <ghps/search_region/local_bounds.hpp>
#include <map>
#include <ghps/mokp/instance.hpp>
#include <ghps/mokp/solvers/branch_bound.hpp>
#include <ghps/mokp/solvers/miqp.hpp>
#include <ghps/mokp/solvers/log.hpp>
#include <mutex>
#include <shared_mutex>
#include <ghps/search_region/rwd.hpp>
#include <ghps/hps/statistics.hpp>
#include <string>
#include <ghps/core/types.hpp>
#include <ghps/search_region/subgraphs.hpp>
#include <ghps/core/timer.hpp>
#include <tuple>
#include <unordered_map>

#include <vector>

namespace moco::ghps {

inline bool DEBUG = false;

enum HPSMode {
  BRANCH_BOUND = 0,
  QUADRATIC_HS = 1,
  LOG_HS = 2,
  HYBRID_HS = 3,
  SIMULATION = 4
};

inline std::vector<Statistics> GHPS(const MOKP& instance, const int32_t K, const int32_t J, const double time_limit = 3600, const int32_t n_threads = 1, HPSMode HPS_MODE = BRANCH_BOUND) {
  if (!instance.is_maximization) throw std::invalid_argument("GHPS optimization supports maximization only");
  if (K <= 0 || K >= std::numeric_limits<int>::digits || J <= 0 || n_threads <= 0
      || !std::isfinite(time_limit) || time_limit < 0 || HPS_MODE < BRANCH_BOUND || HPS_MODE > HYBRID_HS)
    throw std::invalid_argument("Invalid HPS options");
  auto total_time = Timer(time_limit);
  SolverStatus termination = SolverStatus::Optimal;
  std::string termination_message;
  if (DEBUG) {
    std::cout << "Number of objectives: " << instance.M << std::endl;
    std::cout << "Number of items: " << instance.N << std::endl;
    std::cout << "Number of solutions to return: " << J << std::endl;
  }

  std::vector<Solution> solutions;
  std::vector<RefPoint> local_bounds = {RefPoint::origin(instance.M, instance.reference_point)};
  HypervolumeIndicator<hv_t, std::vector<int64_t>> hv_indicator(instance.reference_point, true);

  // Cache for the subproblems - use unordered_map for faster lookups
  std::unordered_map<std::string, MOKPSolution> subproblems_cache;
  std::shared_mutex subproblems_cache_mutex;

  // Statistics
  std::vector<std::tuple<int32_t, int64_t, double>> statistics;
  const std::vector<int32_t> J_values = {5, 10, 15, 20, 50};
  std::vector<Statistics> J_statistics;
  J_statistics.reserve(J_values.size() + 1);

  // Counters
  int64_t n_iterations = 0;
  int64_t n_subproblems = 0;
  int64_t n_cache_hits = 0;
  double subgraphs_time = 0.0;

  while (!local_bounds.empty() && static_cast<int>(solutions.size()) < J) {
    std::cerr << "Iteration " << n_iterations + 1 << ": " << solutions.size() << " solutions found, " << local_bounds.size() << " local bounds remaining." << std::endl;
    // Check if the time limit has been reached
    if (total_time.finished()) {
      termination = SolverStatus::TimeLimit;
      break;
    }
    // Save the statistics for the current size of the solutions if it is in J_values
    if (std::find(J_values.begin(), J_values.end(), static_cast<int>(solutions.size())) != J_values.end()) {
      J_statistics.emplace_back(instance, solutions, statistics, total_time.elapsed(), n_iterations,
                                n_subproblems, n_cache_hits, subgraphs_time, total_time.elapsed() - subgraphs_time);
    }

    if (DEBUG) {
      std::cout << "Iteration: " << n_iterations << std::endl;
      std::cout << "Current bounds:" << std::endl;
      for (const auto& local_bound : local_bounds) {
        std::cout << " - " << local_bound.to_string() << std::endl;
      }
    }

    MOKPSolution global_solution = MOKPSolution();
    std::vector<std::pair<int32_t, int32_t>> graph_edges;
    for (int k = 1; k <= std::min(K, static_cast<int>(local_bounds.size())); k++) {
      // Check if the time limit has been reached
      if (total_time.finished()) {
        break;
      }

      if (DEBUG) {
        std::cout << "Subset size: " << k << std::endl;
      }

      auto subsets_time = Timer();
      std::vector<std::vector<int32_t>> subsets;
      try {
        if (k == 1) {
          subsets = Subgraphs::get_subgraphs(local_bounds, k, total_time);
        } else {
          if (k == 2) graph_edges = Subgraphs::get_graph_edges(local_bounds, total_time);
          subsets = Subgraphs::get_subgraphs(local_bounds.size(), graph_edges, k, total_time);
        }
      } catch (const TimeLimitReached& e) {
        termination = SolverStatus::TimeLimit;
        termination_message = e.what();
        subgraphs_time += subsets_time.elapsed();
        break;
      }
      subgraphs_time += subsets_time.elapsed();

      if (DEBUG) {
        std::cout << "Number of subsets: " << subsets.size() << std::endl;
      }

      MOKPSolution k_solution = MOKPSolution();

      SolverStatus iteration_status = SolverStatus::Optimal;
      std::string iteration_message;
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
        bool cached = false;
        // Protected cache access
        {
          std::shared_lock<std::shared_mutex> lock(subproblems_cache_mutex);
          auto cached_result = subproblems_cache.find(subset_key);
          if (cached_result != subproblems_cache.end()) {
            solution_found = cached_result->second;
            cached = true;
            #pragma omp atomic
            n_cache_hits++;
          }
        }

        if (!cached) {
          int64_t n_recursions = 0;
          auto solution_time = Timer();
          // Prepare data structures for the solvers
          std::vector<std::vector<int64_t>> ref_points_coordinates(ref_points.size());
          for (size_t i = 0; i < ref_points.size(); ++i) {
            ref_points_coordinates[i] = ref_points[i].coordinates;
          }
          MOKPMIQPOptions opt = {.time_limit = std::max(0.0, total_time.left()), .threads = 1, .deadline = total_time};
          try {
            switch (HPS_MODE) {
              case BRANCH_BOUND: {
                auto state = MOKPBranchBound::branch_and_bound_greedy(instance, union_point, ref_points, total_time);
                solution_found = state.solution;
                break;
              }
              case QUADRATIC_HS:
                solution_found = MOKPMIQP(instance, ref_points_coordinates, opt).solve();
                break;
              case LOG_HS:
                solution_found = MOKPLogExp(instance.items, instance.N, instance.M, instance.W,
                    ref_points_coordinates, {.time_limit = std::max(0.0, total_time.left()), .threads = 1, .deadline = total_time}).solve();
                break;
              case HYBRID_HS:
                if (k <= 2)
                  solution_found = MOKPLogExp(instance.items, instance.N, instance.M, instance.W,
                      ref_points_coordinates, {.time_limit = std::max(0.0, total_time.left()), .threads = 1, .deadline = total_time}).solve();
                else
                  solution_found = MOKPMIQP(instance, ref_points_coordinates, opt).solve();
                break;
              default: throw std::invalid_argument("Invalid HPS mode");
            }
          } catch (GRBException& e) {
            solution_found.status = SolverStatus::Error;
            solution_found.message = e.getMessage();
          } catch (const std::exception& e) {
            solution_found.status = SolverStatus::Error;
            solution_found.message = e.what();
          }
          n_recursions = solution_found.nodes;
          solution_found.used_reference_points = ref_points;
          // Second critical section: update cache and stats
          {
            std::unique_lock<std::shared_mutex> lock(subproblems_cache_mutex);
            n_subproblems++;
            if (solution_found.status == SolverStatus::Optimal || solution_found.status == SolverStatus::Infeasible)
              subproblems_cache[subset_key] = solution_found;
            statistics.push_back(std::make_tuple(k, n_recursions, solution_time.elapsed()));
          }
          // std::cout << "Total recursions: " << state.n_recursions << std::endl;
        }

        if (solution_found.status != SolverStatus::Optimal && solution_found.status != SolverStatus::Infeasible) {
          std::lock_guard<std::mutex> lock(k_solution_mutex);
          if (iteration_status == SolverStatus::Optimal || iteration_status == SolverStatus::TimeLimit) {
            iteration_status = solution_found.status;
            iteration_message = std::string(status_name(solution_found.status)) + " (backend status "
                + std::to_string(solution_found.backend_status) + "): " + solution_found.message;
          }
          continue;
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

      if (iteration_status != SolverStatus::Optimal || total_time.finished()) {
        termination = iteration_status == SolverStatus::Optimal ? SolverStatus::TimeLimit : iteration_status;
        termination_message = iteration_message;
        break;
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

    if (termination != SolverStatus::Optimal || total_time.finished()) {
      if (termination == SolverStatus::Optimal) termination = SolverStatus::TimeLimit;
      break;
    }
    if (global_solution.hypervolume <= 0) {
      if (DEBUG) {
        std::cout << "Exiting main loop..." << std::endl;
      }
      break;
    }

    // Compute the real hypervolume of the solution found
    auto real_hypervolume = hv_indicator.insert(global_solution.values);

    Solution solution = Solution("z" + std::to_string(solutions.size()), global_solution.hypervolume, real_hypervolume, global_solution.values);
    solutions.push_back(solution);
    if (DEBUG) {
      std::cout << "Best solution found:" << solution.to_string() << std::endl;
    }

    // Update the local bounds
    try {
      LLB::update_bounds_article(local_bounds, solution, total_time);
    } catch (const TimeLimitReached& e) {
      termination = SolverStatus::TimeLimit;
      termination_message = e.what();
      break;
    }
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

  if (total_time.finished() && termination == SolverStatus::Optimal) termination = SolverStatus::TimeLimit;
  double end_time = total_time.elapsed();
  double subproblems_time = end_time - subgraphs_time;

  // Save the statistics for the last size of the solutions (the J value received)
  J_statistics.emplace_back(instance, solutions, statistics, end_time, n_iterations, n_subproblems, n_cache_hits, subgraphs_time, subproblems_time);

  J_statistics.back().status = termination;
  J_statistics.back().status_message = termination_message;

  if (DEBUG) {
    std::cout << "Greedy HPS finished!" << std::endl;
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

} // namespace moco::ghps

#endif  // GHPS_HPS_ALGORITHM_HPP
