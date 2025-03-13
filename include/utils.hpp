#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <structs.hpp>
#include <vector>
#include <llb.hpp>

class Utils {
 public:
  // Public function to get all subsets of size k from a set of size n
  static std::vector<std::vector<int32_t>> get_k_subsets(int32_t n, int32_t k) {
    std::vector<std::vector<int32_t>> result;
    std::vector<int32_t> combination(k);

    std::function<void(int32_t, int32_t)> backtrack = [&](int32_t start, int32_t index) {
      if (index == k) {
        result.push_back(combination);
        return;
      }
      int end = n - (k - index - 1);
      for (int i = start; i < end; ++i) {
        combination[index] = i;
        backtrack(i + 1, index + 1);
      }
    };
    backtrack(0, 0);
    return result;
  }

  // Public function to compute the exclusion-inclusion principle for hypervolume
  static int64_t compute_hypervolume_exclusion_inclusion(const std::vector<RefPoint>& ref_points, const std::vector<int64_t>& solution) {
    int64_t result = 0;
    int32_t n = ref_points.size();

    for (int32_t i = 1; i < (1 << n); ++i) {
      std::vector<RefPoint> subset;
      for (int32_t j = 0; j < n; ++j) {
        if (i & (1 << j)) {
          subset.push_back(ref_points[j]);
        }
      }

      int64_t hypervolume = compute_hypervolume(subset, solution);
      result += (__builtin_popcount(i) % 2 == 1 ? 1 : -1) * hypervolume;
    }

    return result;
  }

 private:
  // Private function to compute the hypervolume in the branch and bound algorithm
  static int64_t compute_hypervolume(const std::vector<RefPoint>& ref_points, const std::vector<int64_t>& solution) {
    int32_t objectives = ref_points[0].coordinates.size();
    std::vector<int64_t> intersection(objectives, std::numeric_limits<int64_t>::min());
    for (const auto& ref_point : ref_points) {
      for (int32_t m = 0; m < objectives; ++m) {
        intersection[m] = std::max(intersection[m], ref_point.coordinates[m]);
      }
    }
    int64_t hypervolume = 1;
    for (int32_t m = 0; m < objectives; ++m) {
      hypervolume *= solution[m] - intersection[m];
    }
    return hypervolume;
  }
};

#endif  // UTILS_HPP
