#ifndef MOKP_SOLVER_BRANCH_BOUND_HPP
#define MOKP_SOLVER_BRANCH_BOUND_HPP

#include <algorithm>
#include <vector>
#include <iostream>
#include <ghps/core/hypervolume.hpp>
#include <ghps/mokp/instance.hpp>
#include <ghps/core/timer.hpp>

namespace moco::ghps {

struct MOKPBranchBoundState {
  RefPoint union_point;
  std::vector<RefPoint> ref_points;
  std::vector<std::vector<int64_t>> ref_points_coordinates;
  MOKPSolution solution;
  int64_t n_recursions;
  int64_t n_prunings;
  std::vector<int32_t> used_items;

  MOKPBranchBoundState(const RefPoint& up, const std::vector<RefPoint>& rp, int32_t N, int32_t M)
      : union_point(up),
        ref_points(rp),
        solution(0, 0, std::vector<int32_t>(N, 0), std::vector<int64_t>(M, 0)),
        n_recursions(0),
        n_prunings(0),
        used_items(std::vector<int32_t>(N, 0)) {
    // Pre-allocate memory to avoid reallocations
    ref_points_coordinates.reserve(ref_points.size());
    for (const auto& ref_point : ref_points) {
      ref_points_coordinates.push_back(ref_point.coordinates);
    }
  }
};

class MOKPBranchBound {
 private:
  static constexpr bool DEBUG = false;

  static inline bool is_valid(const MOKP& instance, const std::vector<int64_t>& values, const MOKPBranchBoundState& state) {
    // Every coordinate must strictly improve the joint reference.
    for (int i = 0; i < instance.M; i++) {
      if (values[i] <= state.union_point.coordinates[i]) {
        return false;
      }
    }
    return true;
  }

  static void branch_and_bound_greedy_impl(const MOKP& instance, int32_t l, int64_t current_weight,
                               std::vector<int64_t>& values,
                               MOKPBranchBoundState& state,
                               const Timer& timer) {
    if (state.solution.status == SolverStatus::TimeLimit || timer.finished()) {
      state.solution.status = SolverStatus::TimeLimit;
      return;
    }
    state.n_recursions++;
    if (DEBUG) {
      std::cout << "Recursion: " << state.n_recursions << ", l: " << l << ", current_weight: " << current_weight << std::endl;
    }
    // Base case
    if (l == instance.N) {
      if (is_valid(instance, values, state)) {
        hv_t aux_hypervolume = scalarized_hypervolume(state.ref_points_coordinates, values);
        if (aux_hypervolume > state.solution.hypervolume) {
          state.solution = MOKPSolution(aux_hypervolume, current_weight, state.used_items, values);
        }
      }
      return;
    }

    // Compute the upper bound
    std::pair<hv_t, std::vector<int64_t>> upper_bound = compute_bound_improved(instance, l, current_weight, values, state.ref_points_coordinates, timer);
    hv_t upper_bound_hypervolume = upper_bound.first;
    const auto& upper_bound_values = upper_bound.second;
    if (DEBUG) {
      std::cout << "[GREEDY] Greedy upper bound for l = " << l << ": " << upper_bound_hypervolume << ", values: [";
      for (int i = 0; i < instance.M; i++) {
        std::cout << upper_bound_values[i];
        if (i < instance.M - 1) {
          std::cout << ", ";
        }
      }
      std::cout << "]" << std::endl;
    }

    // Prune the branch if the upper bound cannot improve the current solution
    if (upper_bound_hypervolume <= state.solution.hypervolume) {
      state.n_prunings++;
      if (DEBUG) std::cout << "Pruning because upper bound cannot improve the current solution" << std::endl;
      return;
    }
    // Recursive case
    // Include the item if the weight is less than the maximum
    const auto& item = instance.items[instance.sorted_items[0][l]];
    if (item.weight <= instance.W - current_weight) {
      if (DEBUG) std::cout << "Including item " << item.idx << " with weight " << item.weight << std::endl;
      state.used_items[item.idx] = 1;
      timer.check();
      for (int j = 0; j < instance.M; j++) {
        values[j] += item.values[j];
      }
      branch_and_bound_greedy_impl(instance, l + 1, current_weight + item.weight, values, state, timer);
      for (int j = 0; j < instance.M; j++) values[j] -= item.values[j];
      state.used_items[item.idx] = 0;
    }
    if (DEBUG) std::cout << "Excluding item " << item.idx << " with weight " << item.weight << std::endl;
    branch_and_bound_greedy_impl(instance, l + 1, current_weight, values, state, timer);
  }

 public:
  static inline std::pair<hv_t, std::vector<int64_t>> compute_bound_improved(const MOKP& instance, const int32_t& l, const int64_t& current_weight,
                                                                         const std::vector<int64_t>& values,
                                                                         const std::vector<std::vector<int64_t>>& ref_points_coordinates, const Timer& timer = Timer()) {
    std::vector<int64_t> aux_values = values;
    for (int j = 0; j < instance.M; j++) {
      timer.check();
      int64_t aux_weight = current_weight;
      for (int i = 0; i < instance.N; i++) {
        timer.check();
        const auto index = instance.sorted_items[j][i];
        if (instance.item_positions[index] < l) continue;
        const auto& item = instance.items[index];
        if (item.weight <= instance.W - aux_weight) {
          aux_weight += item.weight;
          aux_values[j] += item.values[j];
        } else {
          int64_t remaining_weight = instance.W - aux_weight;
          const hv_t numerator = hv_t(remaining_weight) * item.values[j];
          // Round upward so the fractional bound cannot prune an integer optimum.
          aux_values[j] += static_cast<int64_t>((numerator + item.weight - 1) / item.weight);
          break;
        }
      }
    }
    hv_t aux_hypervolume = scalarized_hypervolume(ref_points_coordinates, aux_values);
    timer.check();
    return std::make_pair(aux_hypervolume, aux_values);
  }

  static Solution simulate_solution(const MOKP& instance, const std::vector<RefPoint>& ref_points, const std::vector<Solution>& curr_sols = {}) {
    auto strongly_dominates = [](const std::vector<int64_t>& v1, const std::vector<int64_t>& v2) {
      for (size_t i = 0; i < v1.size(); ++i) {
        if (v1[i] <= v2[i]) return false;
      }
      return true;
    };

    std::vector<int64_t> union_point = std::vector<int64_t>(instance.M, 0);
    for (size_t i = 0; i < ref_points.size(); ++i) {
      for (int32_t j = 0; j < instance.M; ++j) {
        union_point[j] = std::max(union_point[j], ref_points[i].coordinates[j]);
      }
    }

    std::vector<int64_t> origin = std::vector<int64_t>(instance.M, 0);
    std::vector<std::vector<int64_t>> ref_points_coordinates(ref_points.size());
    for (size_t i = 0; i < ref_points.size(); ++i) {
      ref_points_coordinates[i] = ref_points[i].coordinates;
    }
    Solution best_solution = Solution("z0", 0, 0, origin, ref_points);
    for (const auto& sol : instance.nondominated_set) {
      auto f = std::find_if(curr_sols.begin(), curr_sols.end(),
          [&sol](const Solution& selected) { return selected.coordinates == sol; });
      if (f == curr_sols.end() && strongly_dominates(sol, union_point)) {
        hv_t aux_hypervolume = scalarized_hypervolume(ref_points_coordinates, sol);
        if (aux_hypervolume > best_solution.hypervolume) {
          best_solution.coordinates = sol;
          best_solution.hypervolume = aux_hypervolume;
        }
      }
    }

    return best_solution;
  }

  static Solution compute_greedy_upper_bound(const MOKP& instance, const std::vector<RefPoint>& ref_points) {
    std::vector<std::vector<int64_t>> ref_points_coordinates(ref_points.size());
    for (size_t i = 0; i < ref_points.size(); ++i) {
      ref_points_coordinates[i] = ref_points[i].coordinates;
    }

    auto upper_bound = compute_bound_improved(instance, 0, 0, std::vector<int64_t>(instance.M, 0), ref_points_coordinates);
    hv_t upper_bound_hypervolume = upper_bound.first;
    const auto& upper_bound_values = upper_bound.second;

    Solution best_solution = Solution("z0", upper_bound_hypervolume, 0, upper_bound_values, std::vector<RefPoint>());
    return best_solution;
  }

  static MOKPBranchBoundState branch_and_bound_greedy(const MOKP& instance, const RefPoint& union_point, const std::vector<RefPoint>& ref_points, const Timer& timer) {
    if (!instance.is_maximization) throw std::invalid_argument("Branch-and-bound supports maximization only");
    MOKPBranchBoundState state(union_point, ref_points, instance.N, instance.M);
    std::vector<int64_t> values(instance.M, 0);
    try {
      branch_and_bound_greedy_impl(instance, 0, 0, values, state, timer);
    } catch (const TimeLimitReached& e) {
      state.solution.status = SolverStatus::TimeLimit;
      state.solution.message = e.what();
    }
    if (state.solution.status != SolverStatus::TimeLimit)
      state.solution.status = state.solution.hypervolume > 0 ? SolverStatus::Optimal : SolverStatus::Infeasible;
    state.solution.nodes = state.n_recursions;
    return state;
  }
};

} // namespace moco::ghps

#endif // MOKP_SOLVER_BRANCH_BOUND_HPP
