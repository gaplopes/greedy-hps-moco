#ifndef SIMULATE_HPP
#define SIMULATE_HPP

#include <ghps/core/hypervolume.hpp>
#include <ghps/core/types.hpp>
#include <ghps/core/timer.hpp>
#include <ghps/hps/statistics.hpp>
#include <ghps/mokp/instance.hpp>
#include <ghps/search_region/local_bounds.hpp>

#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace moco::ghps {

template <typename HVType = hv_t>
[[nodiscard]] Statistics GHPS_SIM(const MOKP& instance, int32_t J = -1) {
  if (instance.nondominated_set.empty()) {
    std::cout << "No nondominated solutions found." << std::endl;
    return Statistics(instance, {}, {}, 0.0, 0, 0, 0, 0.0, 0.0);
  }

  const int32_t M = instance.M;
  if (J <= 0 || J > static_cast<int32_t>(instance.nondominated_set.size())) {
    J = static_cast<int32_t>(instance.nondominated_set.size());
  }

  // Pre-convert candidate coordinates once to avoid conversions and narrowing in tight loop
  std::vector<std::vector<HVType>> n_set;
  n_set.reserve(instance.nondominated_set.size());
  for (const auto& point : instance.nondominated_set) n_set.emplace_back(point.begin(), point.end());

  RefPoint ref_point = RefPoint::origin(M, instance.reference_point);
  std::vector<RefPoint> local_bounds = {ref_point};

  std::vector<HVType> ref_coords(ref_point.coordinates.begin(), ref_point.coordinates.end());
  HypervolumeIndicator<HVType, std::vector<HVType>> hvc_space(ref_coords, instance.is_maximization);

  std::vector<Solution> solutions;
  solutions.reserve(static_cast<size_t>(J));

  auto total_time = Timer();

  while (!n_set.empty() && static_cast<int32_t>(solutions.size()) < J) {
    std::cout << "- Current number of solutions: " << solutions.size() << std::endl;
    HVType best_hypervolume = std::numeric_limits<HVType>::lowest();
    int32_t best_solution_index = -1;

    for (size_t i = 0; i < n_set.size(); ++i) {
      HVType sol_hypervolume = hvc_space.contribution(n_set[i]);
      if (sol_hypervolume > 0 && sol_hypervolume > best_hypervolume) {
        best_hypervolume = sol_hypervolume;
        best_solution_index = static_cast<int32_t>(i);
      }
    }

    if (best_solution_index == -1) {
      break;
    }

    // Retrieve best solution and compute dominated bounds only once
    Solution sol("z" + std::to_string(solutions.size()), static_cast<hv_t>(best_hypervolume), 0,
                 std::vector<int64_t>(n_set[best_solution_index].begin(), n_set[best_solution_index].end()));
    sol.dominated_bounds = LLB::get_dominated_bounds(local_bounds, sol);

    hvc_space.insert(n_set[best_solution_index]);
    LLB::update_bounds(local_bounds, sol);
    solutions.push_back(std::move(sol));

    // O(1) swap-and-pop removal from candidates
    std::swap(n_set[best_solution_index], n_set.back());
    n_set.pop_back();
  }

  const double end_time = total_time.elapsed();

  return Statistics(instance, std::move(solutions), {}, end_time, 0, 0, 0, 0.0, 0.0);
}

} // namespace moco::ghps

#endif  // SIMULATE_HPP
