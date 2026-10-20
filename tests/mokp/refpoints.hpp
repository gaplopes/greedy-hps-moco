#ifndef GHPS_EXPERIMENT_REFPOINTS_HPP
#define GHPS_EXPERIMENT_REFPOINTS_HPP

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <ghps/mokp/instance.hpp>
#include <ghps/mokp/solvers/branch_bound.hpp>
#include <ghps/search_region/local_bounds.hpp>

namespace moco::ghps {

inline std::vector<std::vector<int64_t>>
ref_to_v(const std::vector<RefPoint> &ref_points) {
  std::vector<std::vector<int64_t>> ref_points_coordinates(ref_points.size());
  for (size_t i = 0; i < ref_points.size(); ++i) {
    ref_points_coordinates[i] = ref_points[i].coordinates;
  }
  return ref_points_coordinates;
}

inline std::vector<int64_t>
get_union_point(const std::vector<RefPoint> &ref_points) {
  int M = ref_points[0].coordinates.size();
  std::vector<int64_t> union_point(M, 0);
  for (size_t i = 0; i < ref_points.size(); ++i) {
    for (int32_t j = 0; j < M; ++j) {
      union_point[j] = std::max(union_point[j], ref_points[i].coordinates[j]);
    }
  }
  return union_point;
}

inline void generate_reference_points(const std::vector<int64_t>& max_values, int R,
                                      std::vector<RefPoint>& ref_points, int32_t seed,
                                      bool nondominated) {
  if (R <= 0) throw std::invalid_argument("Reference count must be positive");
  std::mt19937 generator(seed);
  std::vector<std::uniform_real_distribution<double>> distributions;
  distributions.reserve(max_values.size());
  for (auto value : max_values) distributions.emplace_back(0.0, value);
  std::set<std::vector<int64_t>> coords_set;
  const int64_t attempt_limit = std::max<int64_t>(1000, int64_t(R) * 1000);
  int64_t attempts = 0;
  while (static_cast<int>(coords_set.size()) < R) {
    if (++attempts > attempt_limit) throw std::runtime_error("Cannot generate requested distinct reference points within attempt limit");
    std::cerr << "Generating reference point " << coords_set.size() + 1 << "/"
              << R << "\r";
    std::vector<int64_t> coordinates;
    coordinates.reserve(max_values.size());
    for (auto& distribution : distributions) coordinates.push_back(distribution(generator));
    // Verify that the generated reference point is not dominated by any
    // reference point in the coordinates set
    if (nondominated && std::any_of(coords_set.begin(), coords_set.end(),
        [&](const auto& existing) { return LLB::weakly_dominates(existing, coordinates); })) continue;
    coords_set.insert(std::move(coordinates));
  }
  ref_points.reserve(ref_points.size() + coords_set.size());
  int rp_idx = 0;
  for (const auto& coords : coords_set) {
    ref_points.emplace_back("r" + std::to_string(rp_idx++), coords);
  }
  std::cerr << "Generated " << ref_points.size()
            << " reference points for testing." << std::endl;
}

inline void generate_reference_points_under_nadir(const MOKP& instance, int R,
                                                 std::vector<RefPoint>& ref_points,
                                                 int32_t& seed) {
  generate_reference_points(instance.get_nadir_point(0), R, ref_points, seed, false);
}

inline void generate_reference_points_anywhere(const MOKP& instance, int R,
                                              std::vector<RefPoint>& ref_points,
                                              int32_t& seed) {
  // Get maximum values from the MOKP::compute_bound_improved
  std::vector<int64_t> values(instance.M, 0);
  RefPoint origin = RefPoint::origin(instance.M);
  const auto max_values =
      MOKPBranchBound::compute_bound_improved(instance, 0, 0, values, {origin.coordinates}).second;
  generate_reference_points(max_values, R, ref_points, seed, false);
}

inline void generate_reference_points_anywhere_nondominated(
    const MOKP& instance, int R, std::vector<RefPoint>& ref_points, int32_t& seed) {
  // Get maximum values from the MOKP::compute_bound_improved
  std::vector<int64_t> values(instance.M, 0);
  RefPoint origin = RefPoint::origin(instance.M);
  const auto max_values =
      MOKPBranchBound::compute_bound_improved(instance, 0, 0, values, {origin.coordinates}).second;
  generate_reference_points(max_values, R, ref_points, seed, true);
}

inline std::vector<std::vector<int32_t>> generate_reference_combinations(
    const std::vector<RefPoint>& ref_points, int max_k, int max_permutations) {
  std::vector<std::vector<int32_t>> combinations;
  for (int k = 2; k <= max_k; k++) {
    // Generate random permutations of k reference points
    std::set<std::vector<int32_t>> tested_permutations;
    std::cerr << "Testing permutations of " << k << " reference points."
              << std::endl;
    const int64_t attempt_limit = std::max<int64_t>(1000, int64_t(max_permutations) * 1000);
    int64_t attempts = 0;
    while (static_cast<int>(tested_permutations.size()) < max_permutations && attempts++ < attempt_limit) {
      std::cerr << "Generating k=" << k << " permutation "
                << tested_permutations.size() + 1 << "/" << max_permutations
                << "\r";
      // Find a set of k reference points that are non-dominated among
      // themselves
      std::vector<int32_t> permutation;
      permutation.reserve(k);
      // Choose one random reference point as the first point
      int32_t first_rp_idx = rand() % ref_points.size();
      permutation.push_back(first_rp_idx);
      for (size_t draw = 0; draw < ref_points.size() * 10 && static_cast<int>(permutation.size()) < k; ++draw) {
        int32_t rp_idx = rand() % ref_points.size();
        // Check if the reference point is non-dominated with respect to the
        // current permutation
        bool dominated = false;
        for (const auto &idx : permutation) {
          if (LLB::weakly_dominates(ref_points[idx].coordinates,
                                    ref_points[rp_idx].coordinates) ||
              LLB::weakly_dominates(ref_points[rp_idx].coordinates,
                                    ref_points[idx].coordinates)) {
            dominated = true;
            break;
          }
        }
        if (!dominated) {
          permutation.push_back(rp_idx);
        }
      }
      if (static_cast<int>(permutation.size()) != k) continue;
      // Test the found permutation
      std::sort(permutation.begin(), permutation.end());
      tested_permutations.insert(std::move(permutation));
    }
    std::cerr << "\n";
    if (static_cast<int>(tested_permutations.size()) < max_permutations)
      throw std::runtime_error("Cannot generate requested incomparable reference combinations for k=" + std::to_string(k)
          + ": generated " + std::to_string(tested_permutations.size()) + " of " + std::to_string(max_permutations));

    combinations.insert(combinations.end(), tested_permutations.begin(), tested_permutations.end());
  }
  return combinations;
}

} // namespace moco::ghps

#endif // GHPS_EXPERIMENT_REFPOINTS_HPP
