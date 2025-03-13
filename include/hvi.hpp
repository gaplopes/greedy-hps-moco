#ifndef HVI_HPP
#define HVI_HPP

#include <cstdint>
#include <vector>

#include <mooutils/indicators.hpp>

template <typename T, typename Solution>
class HypervolumeIndicator {
 public:
  HypervolumeIndicator(Solution ref_point, bool is_maximization)
      : ref_point(ref_point),
        is_maximization(is_maximization),
        hv_space(initialize_hv_space(ref_point, is_maximization)) {}

  T insert(const Solution& sol) {
    Solution aux_sol = is_maximization ? sol : negative_vector(sol);
    return hv_space.insert(aux_sol);
  }

  T contribution(const Solution& sol) const {
    Solution aux_sol = is_maximization ? sol : negative_vector(sol);
    return hv_space.contribution(aux_sol);
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

  T set_contribution(const std::vector<Solution>& solutions) {
    T hvc = 0;
    for (const Solution& sol : solutions) {
      hvc += contribution(sol);
    }
    return hvc;
  }

 private:
  static mooutils::incremental_hv<T, Solution> initialize_hv_space(const Solution& ref_point, bool is_maximization) {
    return is_maximization ? mooutils::incremental_hv<T, Solution>(ref_point)
                           : mooutils::incremental_hv<T, Solution>(negative_vector(ref_point));
  }

  static Solution negative_vector(const Solution& v) {
    Solution res(v.size());
    std::transform(v.begin(), v.end(), res.begin(), std::negate<T>());
    return res;
  }

  const Solution ref_point;
  const bool is_maximization;
  mooutils::incremental_hv<T, Solution> hv_space;
};

#endif  // HVI_HPP
