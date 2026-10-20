#ifndef GHPS_CORE_HYPERVOLUME_HPP
#define GHPS_CORE_HYPERVOLUME_HPP

#include <cstdint>
#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include <mooutils/indicators.hpp>
#include <ghps/core/types.hpp>

namespace moco::ghps {


template <typename T, typename Solution>
class HypervolumeIndicator {
 public:
  HypervolumeIndicator(const Solution& ref_point, bool is_maximization)
      : is_maximization(is_maximization),
        reference_point(ref_point),
        hv_space(initialize_hv_space(ref_point, is_maximization)) {}

  T insert(const Solution& sol) {
    if (!dominates_reference(sol)) return 0;
    return is_maximization ? hv_space.insert(sol) : hv_space.insert(negative_vector(sol));
  }

  T contribution(const Solution& sol) const {
    if (!dominates_reference(sol)) return 0;
    return is_maximization ? hv_space.contribution(sol) : hv_space.contribution(negative_vector(sol));
  }

  T value() const {
    return hv_space.value();
  }

  T set_hypervolume(const std::vector<Solution>& solutions) {
    T hvc = 0;
    for (const Solution& sol : solutions) {
      hvc += insert(sol);
    }
    return hvc;
  }

  T set_hypervolume_unchecked(const std::vector<Solution>& solutions) {
    return set_hypervolume(solutions);
  }

  static T compute_contribution(const std::vector<Solution>& ref_points, const Solution& sol) {
    // Hypervolume of the region bounded below by ref_points and above by sol
    // Using minimization mode (false) negates the coordinates, transforming the problem
    // into finding the volume dominated by -ref_points bounded below by -sol.
    return HypervolumeIndicator<T, Solution>(sol, false).set_hypervolume(ref_points);
  }

 private:
  bool dominates_reference(const Solution& sol) const {
    if (sol.size() != reference_point.size()) throw std::invalid_argument("Hypervolume dimension mismatch");
    for (size_t i = 0; i < sol.size(); ++i) {
      if (is_maximization ? sol[i] <= reference_point[i] : sol[i] >= reference_point[i]) return false;
    }
    return true;
  }

  static mooutils::incremental_hv<T, Solution> initialize_hv_space(const Solution& ref_point, bool is_maximization) {
    if (ref_point.size() < 2) throw std::invalid_argument("Hypervolume requires at least two objectives");
    return is_maximization ? mooutils::incremental_hv<T, Solution>(ref_point)
                           : mooutils::incremental_hv<T, Solution>(negative_vector(ref_point));
  }

  static Solution negative_vector(const Solution& v) {
    Solution res(v.size());
    std::transform(v.begin(), v.end(), res.begin(), std::negate<T>());
    return res;
  }

  const bool is_maximization;
  const Solution reference_point;
  mooutils::incremental_hv<T, Solution> hv_space;
};

inline hv_t scalarized_hypervolume(const std::vector<std::vector<int64_t>>& refs,
                                   const std::vector<int64_t>& values) {
  if (refs.empty() || values.empty()) return 0;
  std::vector<hv_t> lower(values.begin(), values.end());
  std::vector<std::vector<hv_t>> wide_refs;
  for (const auto& ref : refs) {
    if (ref.size() != values.size()) throw std::invalid_argument("Reference dimension mismatch");
    for (size_t d = 0; d < values.size(); ++d) {
      if (values[d] <= ref[d]) return 0;
      lower[d] = std::min(lower[d], hv_t(ref[d]));
    }
    wide_refs.emplace_back(ref.begin(), ref.end());
  }
  hv_t box = 1;
  for (size_t d = 0; d < values.size(); ++d) {
    const hv_t extent = hv_t(values[d]) - lower[d];
    if (box > std::numeric_limits<hv_t>::max() / extent)
      throw std::overflow_error("Scalarized hypervolume exceeds hv_t range");
    box *= extent;
  }
  // Strict dominance was checked above; redundant references still form a union.
  return HypervolumeIndicator<hv_t, std::vector<hv_t>>(std::vector<hv_t>(values.begin(), values.end()), false)
      .set_hypervolume_unchecked(wide_refs);
}

} // namespace moco::ghps

#endif  // GHPS_CORE_HYPERVOLUME_HPP
