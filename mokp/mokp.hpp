#ifndef MOKP_HPP
#define MOKP_HPP

#include <algorithm>
#include <cstdint>
#include <hvi.hpp>
#include <iostream>
#include <structs.hpp>
#include <time.hpp>
#include <utils.hpp>
#include <vector>

struct Item {
 public:
  int32_t idx;                  // Index of the item
  int64_t weight;               // Weight of the item
  std::vector<int64_t> values;  // Values of the item for each objective

  Item(int idx, int64_t weight, std::vector<int64_t> values)
      : idx(idx), weight(weight), values(values) {}

  std::string to_string() const {
    std::string str = "Item(" + std::to_string(idx) + ", " + std::to_string(weight) + ", [";
    for (int i = 0; i < (int)values.size(); i++) {
      str += std::to_string(values[i]);
      if (i < (int)values.size() - 1) {
        str += ", ";
      }
    }
    str += "])\n";
    return str;
  }
};

struct MOKPSolution {
 public:
  int64_t hypervolume;
  int64_t weight;
  std::vector<int32_t> items;
  std::vector<int64_t> values;
  std::vector<RefPoint> used_reference_points;

  MOKPSolution(int64_t hypervolume = 0, int64_t weight = 0, const std::vector<int32_t>& items = {}, const std::vector<int64_t>& values = {}, const std::vector<RefPoint>& used_reference_points = {})
      : hypervolume(hypervolume), weight(weight), items(items), values(values), used_reference_points(used_reference_points) {}

  std::string to_string() const {
    std::string str = "(" + std::to_string(hypervolume) + ", " + std::to_string(weight) + ", [";
    for (int i = 0; i < (int)items.size(); i++) {
      str += std::to_string(items[i]);
      if (i < (int)items.size() - 1) {
        str += ", ";
      }
    }
    str += "], [";
    for (int i = 0; i < (int)values.size(); i++) {
      str += std::to_string(values[i]);
      if (i < (int)values.size() - 1) {
        str += ", ";
      }
    }
    str += "])\n";
    str += "Used reference points: [";
    for (int i = 0; i < (int)used_reference_points.size(); i++) {
      str += used_reference_points[i].to_string() + ", "[i < (int)used_reference_points.size() - 1];
    }
    str += "]\n";
    return str;
  };
};

struct MOKPRecursiveState {
  RefPoint union_point;
  std::vector<RefPoint> ref_points;
  std::vector<std::vector<int64_t>> ref_points_coordinates;
  MOKPSolution solution;
  int64_t n_recursions;
  std::vector<int32_t> used_items;

  MOKPRecursiveState(const RefPoint& up, const std::vector<RefPoint>& rp, int32_t N, int32_t M)
      : union_point(up),
        ref_points(rp),
        solution(0, 0, std::vector<int32_t>(N, 0), std::vector<int64_t>(M, 0)),
        n_recursions(0),
        used_items(std::vector<int32_t>(N, 0)) {
    ref_points_coordinates.reserve(ref_points.size());
    for (const auto& ref_point : ref_points) {
      ref_points_coordinates.push_back(ref_point.coordinates);
    }
  }
};

struct MOKP {
  int32_t N;                                    // Number of elements for each objective
  int32_t M;                                    // Number of objectives
  int64_t W;                                    // Maximum Knapsack weight
  std::vector<int64_t> weights;                 // Weights of each element
  std::vector<std::vector<int64_t>> values;     // Values of each element for each objective
  std::vector<Item> items;                      // Items
  std::vector<std::vector<Item>> sorted_items;  // Sorted items

  std::vector<int64_t> reference_point;    // Reference point
  int32_t n_nondominated_solutions;        // Total number of solutions in the nondominated set
  std::vector<Solution> nondominated_set;  // Nondominated set

  MOKP(int32_t N, int32_t M, int64_t W,
       std::vector<int64_t> weights,
       std::vector<std::vector<int64_t>> values,
       std::vector<Item> items,
       std::vector<Solution> nondominated_set)
      : N(N),
        M(M),
        W(W),
        weights(weights),
        values(values),
        items(items),
        n_nondominated_solutions(nondominated_set.size()),
        nondominated_set(nondominated_set) {
    // Initialize the reference point
    this->reference_point = std::vector<int64_t>(M, 0);
    // Sort the items by the ratio of values[i]/weight lexicographically
    this->sorted_items.resize(M);
    for (int i = 0; i < M; i++) {
      this->sorted_items[i] = (i == 0) ? items : this->sorted_items[0];
      // Sort the items by the ratio of values[i]/weight
      std::sort(sorted_items[i].begin(), sorted_items[i].end(), [i](const Item& a, const Item& b) {
        return (a.values[i] / (double)a.weight) > (b.values[i] / (double)b.weight);
      });
      // If i == 0, rewrite the idx of the items
      if (i == 0) {
        for (int j = 0; j < N; j++) {
          sorted_items[i][j].idx = j;
        }
      }
    }
  }

  template <typename IStream>
  static auto from_stream(IStream&& is) -> MOKP {
    // Read from file
    int32_t N, M;
    int64_t W;
    is >> N >> M;
    is >> W;
    std::vector<int64_t> weights(N);
    std::vector<std::vector<int64_t>> values(N, std::vector<int64_t>(M));
    for (int i = 0; i < N; i++) {
      is >> weights[i];
      for (int j = 0; j < M; j++) {
        is >> values[i][j];
      }
    }
    // Pre-process the items and sort them by the ratio of values[i]/weight
    std::vector<Item> items = std::vector<Item>();
    for (int i = 0; i < N; i++) {
      items.push_back(Item(i, weights[i], values[i]));
    }
    // Read the nondominated set
    std::vector<Solution> nondominated_set;
    int32_t n_nondominated_set = 0;
    // Try to read the number of nondominated solutions
    if (is >> n_nondominated_set) {
      nondominated_set = std::vector<Solution>(n_nondominated_set);
      for (int i = 0; i < n_nondominated_set; i++) {
        std::vector<int64_t> coordinates(M);
        for (int j = 0; j < M; j++) {
          is >> coordinates[j];
        }
        nondominated_set[i] = Solution("z" + std::to_string(i), 0, 0, coordinates, std::vector<RefPoint>());
      }
    }
    std::cout << "MOKP instance loaded successfully!" << std::endl;
    return MOKP(N, M, W, weights, values, items, nondominated_set);
  }

  MOKPRecursiveState compute_solution(const RefPoint& union_point, const std::vector<RefPoint>& ref_points, const Timer& timer) {
    MOKPRecursiveState state(union_point, ref_points, N, M);
    compute_upper_bound_solution(0, 0, std::vector<int64_t>(M, 0), state, timer);
    return state;
  }

 private:
  void compute_upper_bound_solution(int32_t l, int64_t current_weight,
                                    std::vector<int64_t> values,
                                    MOKPRecursiveState& state,
                                    const Timer& timer) {
    if (timer.finished()) {
      return;
    }
    state.n_recursions++;
    // Base case
    if (is_valid(l, values, state)) {
      // int64_t aux_hypervolume = Utils::compute_hypervolume_exclusion_inclusion(state.ref_points, values);
      int64_t aux_hypervolume = HypervolumeIndicator<int64_t, std::vector<int64_t>>(values, false).set_hypervolume(state.ref_points_coordinates);
      if (aux_hypervolume > state.solution.hypervolume) {
        state.solution = MOKPSolution(aux_hypervolume, current_weight, state.used_items, values);
      }
    }
    // No more items to consider
    if (l == N) {
      return;
    }
    // Compute the upper bound
    std::pair<int64_t, std::vector<int64_t>> upper_bound = compute_bound_improved(l, current_weight, values, state);
    int64_t upper_bound_hypervolume = upper_bound.first;
    std::vector<int64_t> upper_bound_values = upper_bound.second;
    // Prune the branch if the upper bound is worse than the current solution
    if (upper_bound_hypervolume < state.solution.hypervolume) {
      // std::cout << "Pruning because upper bound is worse than the current solution" << std::endl;
      return;
    }
    // Prune the branch if the upper bound values are worse than the union point
    if (is_valid(l, upper_bound_values, state) == false) {
      // std::cout << "Pruning because upper bound values are worse than the union point" << std::endl;
      return;
    }
    // Recursive case
    int aux_weight = sorted_items[0][l].weight + current_weight;
    // Include the item if the weight is less than the maximum
    if (aux_weight <= this->W) {
      state.used_items[l] = 1;
      std::vector<int64_t> aux_values = values;
      for (int j = 0; j < this->M; j++) {
        aux_values[j] += sorted_items[0][l].values[j];
      }
      compute_upper_bound_solution(l + 1, current_weight + sorted_items[0][l].weight, aux_values, state, timer);
      state.used_items[l] = 0;
    }
    // Exclude the item
    compute_upper_bound_solution(l + 1, current_weight, values, state, timer);
  }

  std::pair<int64_t, std::vector<int64_t>> compute_bound_improved(const int32_t& l, const int64_t& current_weight,
                                                                  const std::vector<int64_t>& values,
                                                                  const MOKPRecursiveState& state) {
    std::vector<int64_t> aux_values = values;
    for (int j = 0; j < this->M; j++) {
      int64_t aux_weight = current_weight;
      for (int i = 0; i < N; i++) {
        if (j == 0 && i < l) continue;
        if (j != 0 && sorted_items[j][i].idx < l) continue;
        if (aux_weight + sorted_items[j][i].weight <= this->W) {
          aux_weight += sorted_items[j][i].weight;
          aux_values[j] += sorted_items[j][i].values[j];
        } else {
          int64_t remaining_weight = this->W - aux_weight;
          double fraction = remaining_weight / (double)sorted_items[j][i].weight;
          aux_values[j] += fraction * sorted_items[j][i].values[j];
          aux_weight += fraction * sorted_items[j][i].weight;
          break;
        }
      }
    }
    // int64_t aux_hypervolume = Utils::compute_hypervolume_exclusion_inclusion(state.ref_points, aux_values);
    int64_t aux_hypervolume = HypervolumeIndicator<int64_t, std::vector<int64_t>>(aux_values, false).set_hypervolume(state.ref_points_coordinates);
    return std::make_pair(aux_hypervolume, aux_values);
  }

  inline bool is_valid(const int32_t& l, const std::vector<int64_t>& values, const MOKPRecursiveState& state) {
    bool better = true;
    better = better && l <= N;
    for (int i = 0; i < this->M; i++) {
      better = better && (values[i] > state.union_point.coordinates[i]);
    }
    return better;
  }
};

#endif  // MOKP_HPP
