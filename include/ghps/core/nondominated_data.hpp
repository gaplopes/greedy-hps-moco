#ifndef GHPS_CORE_NONDOMINATED_DATA_HPP
#define GHPS_CORE_NONDOMINATED_DATA_HPP

#include <algorithm>
#include <cstdint>
#include <ghps/core/hypervolume.hpp>
#include <ghps/core/types.hpp>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace moco::ghps {

namespace detail {

struct NondominatedData {
  const bool is_maximization;
  const std::vector<std::vector<int64_t>> points;
  const std::set<std::vector<int64_t>> membership;
  const std::vector<int64_t> reference;
  const std::vector<int64_t> nadir;
  const hv_t hypervolume;
  const hv_t nadir_hypervolume;

  NondominatedData(int32_t dimensions, std::vector<std::vector<int64_t>> points, bool is_maximization)
      : is_maximization(is_maximization), points(std::move(points)), membership(this->points.begin(), this->points.end()),
        reference(dimensions, 0), nadir(compute_nadir()),
        hypervolume(compute_hypervolume(reference)), nadir_hypervolume(compute_hypervolume(nadir)) {}

 private:
  std::vector<int64_t> compute_nadir() const {
    if (points.empty()) return reference;
    auto result = points.front();
    for (const auto& point : points)
      for (size_t i = 0; i < result.size(); ++i)
        result[i] = is_maximization ? std::min(result[i], point[i]) : std::max(result[i], point[i]);
    for (auto& value : result) {
      if (!is_maximization && value == std::numeric_limits<int64_t>::max())
        throw std::overflow_error("Nadir coordinate exceeds int64_t range");
      value += is_maximization ? -1 : 1;
    }
    return result;
  }

  hv_t compute_hypervolume(const std::vector<int64_t>& reference) const {
    if (points.empty()) return 0;
    std::vector<std::vector<hv_t>> coordinates;
    coordinates.reserve(points.size());
    for (const auto& point : points) coordinates.emplace_back(point.begin(), point.end());
    return HypervolumeIndicator<hv_t, std::vector<hv_t>>(
               std::vector<hv_t>(reference.begin(), reference.end()), is_maximization)
        .set_hypervolume_unchecked(coordinates);
  }
};

} // namespace detail

} // namespace moco::ghps

#endif // GHPS_CORE_NONDOMINATED_DATA_HPP
