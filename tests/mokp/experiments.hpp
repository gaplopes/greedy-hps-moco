#ifndef GHPS_EXPERIMENTS_HPP
#define GHPS_EXPERIMENTS_HPP

#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "config.hpp"
#include "refpoints.hpp"
#include "solver.hpp"

namespace moco::ghps {

inline void experiments(const MOKP& instance, const TestConfig& config) {
  if (config.num_ref_points <= 0 || config.max_k <= 0 || config.max_k > config.num_ref_points
      || config.max_k >= std::numeric_limits<int>::digits || config.max_permutations <= 0
      || config.solver_threads <= 0)
    throw std::invalid_argument("Invalid experiment options: require positive counts, max-k <= refs and max-k < 31");
  if (instance.nondominated_set.empty()) throw std::invalid_argument("Experiments require a supplied nondominated set");
  srand(config.seed);
  ExperimentRunner runner(instance, config);

  // 1 - Generate random reference points
  int R = config.num_ref_points; // Number of reference points
  std::vector<RefPoint> ref_points;
  int32_t seed = config.seed;
  // generate_reference_points_under_nadir(instance, R, ref_points, seed);
  generate_reference_points_anywhere(instance, R, ref_points, seed);
  // Add origin as a reference point if not already present
  // ref_points.push_back(RefPoint::origin(instance.M));

  const auto combinations = instance.M == 2 ? std::vector<std::vector<int32_t>>{}
      : generate_reference_combinations(ref_points, config.max_k, config.max_permutations);

  // 2.1 - Individual reference points
  for (const auto& rp : ref_points) runner.run({rp});

  // 2.2 - Combinations of reference points
  for (const auto& combination : combinations) {
    std::vector<RefPoint> current_ref_points;
    current_ref_points.reserve(combination.size());
    for (auto index : combination) current_ref_points.push_back(ref_points[index]);
    runner.run(current_ref_points);
  }
  runner.finish();
  if (instance.M == 2)
    std::cerr << "Skipping combinations of reference points for M=2." << std::endl;
  std::cerr << "Sequential execution finished." << std::endl;
}

} // namespace moco::ghps

#endif // GHPS_EXPERIMENTS_HPP
