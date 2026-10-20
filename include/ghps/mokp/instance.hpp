#ifndef MOKP_INSTANCE_HPP
#define MOKP_INSTANCE_HPP

#include <algorithm>
#include <cstdint>
#include <ghps/core/nondominated_data.hpp>
#include <ghps/core/types.hpp>
#include <ghps/mokp/item.hpp>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace moco::ghps {

struct MOKPSolution {
 public:
  hv_t hypervolume;
  int64_t weight;
  std::vector<int32_t> items;  // Binary selections in original input row order
  std::vector<int64_t> values;
  std::vector<RefPoint> used_reference_points;
  SolverStatus status = SolverStatus::NotSolved;
  int backend_status = 0;
  int64_t nodes = 0;
  std::string message;

  MOKPSolution(hv_t hypervolume = 0,
               int64_t weight = 0,
               std::vector<int32_t> items = {},
               std::vector<int64_t> values = {},
               std::vector<RefPoint> used_reference_points = {})
      : hypervolume(hypervolume),
        weight(weight),
        items(std::move(items)),
        values(std::move(values)),
        used_reference_points(std::move(used_reference_points)) {}

  std::string to_string() const {
    std::stringstream ss;
    ss << "(" << hypervolume << ", " << weight << ", [";
    for (int i = 0; i < (int)items.size(); i++) {
      ss << items[i];
      if (i < (int)items.size() - 1) {
        ss << ", ";
      }
    }
    ss << "], [";
    for (int i = 0; i < (int)values.size(); i++) {
      ss << values[i];
      if (i < (int)values.size() - 1) {
        ss << ", ";
      }
    }
    ss << "])\n";
    ss << "Used reference points: [";
    for (int i = 0; i < (int)used_reference_points.size(); i++) {
      ss << used_reference_points[i].to_string();
      if (i < (int)used_reference_points.size() - 1) {
        ss << ", ";
      }
    }
    ss << "]\n";
    ss << "Status: " << status_name(status) << ", backend status: " << backend_status << "\n";
    if (!message.empty()) ss << message << "\n";
    return ss.str();
  };
};

struct MOKP {
  const int32_t N;                                           // Number of items
  const int32_t M;                                           // Number of objectives
  const int64_t W;                                           // Maximum Knapsack weight
  const std::vector<Item> items;                             // Items in original input row order
  const std::vector<std::vector<int32_t>> sorted_items;      // Input row indices sorted by each objective ratio
  const std::vector<int32_t> item_positions;                 // Search position of each original input row

 private:
  const std::shared_ptr<const detail::NondominatedData> front_data;

 public:
  const bool& is_maximization;
  const std::vector<int64_t>& reference_point;               // Reference point
  const int32_t n_nondominated_solutions;                    // Total number of solutions in the nondominated set
  const std::vector<std::vector<int64_t>>& nondominated_set; // Nondominated set

  const std::set<std::vector<int64_t>>& nd_set_values;

  MOKP(int32_t N, int32_t M, int64_t W,
       std::vector<int64_t> weights,
       std::vector<std::vector<int64_t>> values,
       std::vector<Item> items,
       std::vector<std::vector<int64_t>> nondominated_set, bool is_maximization = true)
      : MOKP(N, M, W, checked_items(weights, values, std::move(items)), std::move(nondominated_set), is_maximization) {}

  MOKP(int32_t N, int32_t M, int64_t W,
       std::vector<Item> items,
       std::vector<std::vector<int64_t>> nondominated_set, bool is_maximization = true)
      : N(validate(N, M, W, items, nondominated_set)),
        M(M),
        W(W),
        items(input_order(std::move(items))),
        sorted_items(sort_items(this->items, M)),
        item_positions(compute_item_positions(sorted_items.front())),
        front_data(std::make_shared<const detail::NondominatedData>(M, std::move(nondominated_set), is_maximization)),
        is_maximization(front_data->is_maximization),
        reference_point(front_data->reference),
        n_nondominated_solutions(front_data->points.size()),
        nondominated_set(front_data->points),
        nd_set_values(front_data->membership),
        upper_bound(compute_dantzig_upper_bound()) {}
  
  hv_t get_nondominated_set_hv() const { return front_data->hypervolume; }
  hv_t get_nadir_nset_hv() const { return front_data->nadir_hypervolume; }
  std::shared_ptr<const detail::NondominatedData> reporting_data() const { return front_data; }

  /// Dantzig (fractional knapsack) upper bound vector, computed at construction.
  const std::vector<double>& dantzig_upper_bound() const { return upper_bound; }

  std::vector<int64_t> get_nadir_point(int64_t epsilon = 1) const {
    auto result = front_data->nadir;
    if (!nondominated_set.empty()) {
      for (auto& value : result) {
        const hv_t adjusted = hv_t(value) + (is_maximization ? 1 - hv_t(epsilon) : hv_t(epsilon) - 1);
        if (adjusted < std::numeric_limits<int64_t>::min() || adjusted > std::numeric_limits<int64_t>::max())
          throw std::overflow_error("Nadir coordinate exceeds int64_t range");
        value = static_cast<int64_t>(adjusted);
      }
    }
    return result;
  }

 private:
  const std::vector<double> upper_bound;

  static int32_t validate(int32_t N, int32_t M, int64_t W,
                          const std::vector<Item>& items,
                          const std::vector<std::vector<int64_t>>& nondominated_set) {
    if (N <= 0 || M < 2 || M >= std::numeric_limits<int>::digits || W < 0)
      throw std::invalid_argument("Expected N > 0, 2 <= M < 31, and nonnegative capacity");
    if (items.size() != size_t(N))
      throw std::invalid_argument("Item count does not match N");
    std::vector<hv_t> totals(M, 0);
    std::vector<bool> seen(N, false);
    for (int i = 0; i < N; ++i) {
      if (items[i].idx < 0 || items[i].idx >= N || seen[items[i].idx])
        throw std::invalid_argument("Item IDs must be a permutation of the input row indices");
      seen[items[i].idx] = true;
      if (items[i].weight < 0 || items[i].values.size() != size_t(M))
        throw std::invalid_argument("Invalid item weight or objective count");
      for (int j = 0; j < M; ++j) {
        if (items[i].values[j] < 0)
          throw std::invalid_argument("Objective values must be nonnegative");
        totals[j] += items[i].values[j];
        if (totals[j] > std::numeric_limits<int64_t>::max())
          throw std::invalid_argument("Objective sum exceeds the coordinate range");
      }
    }
    for (const auto& point : nondominated_set) {
      if (point.size() != size_t(M) || std::any_of(point.begin(), point.end(), [](int64_t v) { return v < 0; }))
        throw std::invalid_argument("Invalid nondominated coordinates");
    }
    return N;
  }

  static std::vector<Item> checked_items(const std::vector<int64_t>& weights,
                                        const std::vector<std::vector<int64_t>>& values,
                                        std::vector<Item> items) {
    if (weights.size() != items.size() || values.size() != items.size())
      throw std::invalid_argument("Item count does not match weight and objective rows");
    for (const auto& item : items) {
      if (item.idx < 0 || size_t(item.idx) >= items.size()
          || item.weight != weights[item.idx] || item.values != values[item.idx])
        throw std::invalid_argument("Items disagree with input weight or objective rows");
    }
    return items;
  }

  static std::vector<Item> input_order(std::vector<Item> items) {
    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.idx < b.idx; });
    return items;
  }

  static std::vector<std::vector<int32_t>> sort_items(const std::vector<Item>& items, int32_t M) {
    std::vector<std::vector<int32_t>> sorted(M);
    for (int i = 0; i < M; i++) {
      if (i == 0) {
        sorted[i].resize(items.size());
        std::iota(sorted[i].begin(), sorted[i].end(), 0);
      } else {
        sorted[i] = sorted[0];
      }
      // Sort the items by the ratio of values[i]/weight
      std::sort(sorted[i].begin(), sorted[i].end(), [&items, i](int32_t left, int32_t right) {
        const auto& a = items[left];
        const auto& b = items[right];
        if (a.weight == 0 || b.weight == 0) {
          if (a.weight == 0 && b.weight == 0) return (a.values[i] > 0) > (b.values[i] > 0);
          if (a.weight == 0) return a.values[i] > 0;
          return b.values[i] == 0 && a.values[i] > 0;
        }
        return hv_t(a.values[i]) * b.weight > hv_t(b.values[i]) * a.weight;
      });
    }
    return sorted;
  }

  static std::vector<int32_t> compute_item_positions(const std::vector<int32_t>& order) {
    std::vector<int32_t> positions(order.size());
    for (size_t i = 0; i < order.size(); ++i) positions[order[i]] = i;
    return positions;
  }

  std::vector<double> compute_dantzig_upper_bound() const {
    std::vector<double> result(M);
    for (int j = 0; j < M; ++j) {
      int64_t remaining = W;
      for (const auto index : sorted_items[j]) {
        const auto& item = items[index];
        if (remaining <= 0 && item.weight != 0) break;
        if (item.weight <= remaining) {
          result[j] += item.values[j];
          remaining -= item.weight;
        } else {
          result[j] += item.values[j] * (static_cast<double>(remaining) / item.weight);
          remaining = 0;
        }
      }
    }
    return result;
  }

 public:
  template <typename IStream>
  static auto from_stream(IStream&& is, bool is_maximization = true) -> MOKP {
    // Read from file
    int32_t N, M;
    int64_t W;
    if (!(is >> N >> M >> W) || N <= 0 || M < 2 || M >= std::numeric_limits<int>::digits || W < 0)
      throw std::invalid_argument("Invalid instance header: expected N > 0, 2 <= M < 31, W >= 0");
    std::vector<Item> items;
    items.reserve(N);
    for (int i = 0; i < N; i++) {
      int64_t weight;
      if (!(is >> weight) || weight < 0) throw std::invalid_argument("Invalid or missing item weight");
      std::vector<int64_t> row(M);
      for (int j = 0; j < M; j++) {
        if (!(is >> row[j]) || row[j] < 0) throw std::invalid_argument("Invalid or missing item objective");
      }
      items.emplace_back(i, weight, std::move(row));
    }
    // Read the nondominated set
    std::vector<std::vector<int64_t>> nondominated_set;
    int32_t n_nondominated_set = 0;
    // Try to read the number of nondominated solutions
    is >> std::ws;
    if (!is.eof()) {
      if (!(is >> n_nondominated_set) || n_nondominated_set < 0)
        throw std::invalid_argument("Invalid nondominated solution count");
      for (int i = 0; i < n_nondominated_set; i++) {
        std::vector<int64_t> coordinates(M);
        for (int j = 0; j < M; j++) {
          if (!(is >> coordinates[j]) || coordinates[j] < 0)
            throw std::invalid_argument("Invalid or missing nondominated coordinate");
        }
        nondominated_set.push_back(std::move(coordinates));
      }
      is >> std::ws;
      if (!is.eof()) throw std::invalid_argument("Unexpected data after nondominated set");
    }
    MOKP instance(N, M, W, std::move(items), std::move(nondominated_set), is_maximization);
    std::cerr << "MOKP instance loaded successfully!" << std::endl;
    return instance;
  }
};

}  // namespace moco::ghps

#endif  // MOKP_INSTANCE_HPP
