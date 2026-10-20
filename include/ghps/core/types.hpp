#ifndef GHPS_CORE_TYPES_HPP
#define GHPS_CORE_TYPES_HPP

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace moco::ghps {

__extension__ typedef __int128 int128_t;
__extension__ typedef unsigned __int128 uint128_t;
using hv_t = int128_t;

enum class SolverStatus { NotSolved, Optimal, Infeasible, TimeLimit, Interrupted, NumericalFailure, Unbounded, Error, Other };

inline const char* status_name(SolverStatus status) {
  switch (status) {
    case SolverStatus::NotSolved: return "not_solved";
    case SolverStatus::Optimal: return "optimal";
    case SolverStatus::Infeasible: return "infeasible";
    case SolverStatus::TimeLimit: return "time_limit";
    case SolverStatus::Interrupted: return "interrupted";
    case SolverStatus::NumericalFailure: return "numerical_failure";
    case SolverStatus::Unbounded: return "unbounded";
    case SolverStatus::Error: return "error";
    default: return "other";
  }
}

inline std::string to_string(int128_t val) {
  if (val == 0) return "0";
  bool neg = val < 0;
  if (neg) val = -val;
  std::string s;
  while (val > 0) {
    s += static_cast<char>('0' + (val % 10));
    val /= 10;
  }
  if (neg) s += '-';
  std::reverse(s.begin(), s.end());
  return s;
}

inline std::string to_string(uint128_t val) {
  if (val == 0) return "0";
  std::string s;
  while (val > 0) {
    s += static_cast<char>('0' + (val % 10));
    val /= 10;
  }
  std::reverse(s.begin(), s.end());
  return s;
}

inline std::ostream& operator<<(std::ostream& os, int128_t val) {
  return os << to_string(val);
}

inline std::ostream& operator<<(std::ostream& os, uint128_t val) {
  return os << to_string(val);
}

struct RefPoint {
 public:
  std::string id;
  std::vector<int64_t> coordinates;
  std::vector<std::string> defining_points;

  RefPoint() = default;

  RefPoint(std::string id, std::vector<int64_t> coordinates, std::vector<std::string> defining_points = {})
      : id(std::move(id)), coordinates(std::move(coordinates)), defining_points(std::move(defining_points)) {}

  bool operator==(const RefPoint& other) const {
    return coordinates == other.coordinates;
  }

  static RefPoint origin(int32_t M, std::vector<int64_t> coordinates = {}) {
    if (coordinates.empty()) {
      coordinates = std::vector<int64_t>(M, 0);
    }
    std::string id = "u0";
    std::vector<std::string> defining_points;
    defining_points.reserve(M);
    for (int i = 0; i < M; i++) {
      defining_points.push_back("d" + std::to_string(i + 1));
    }
    return RefPoint(std::move(id), std::move(coordinates), std::move(defining_points));
  }

  std::string to_string() const {
    std::string s = id + " (";
    for (size_t i = 0; i < coordinates.size(); ++i) {
      s += std::to_string(coordinates[i]);
      if (i != coordinates.size() - 1)
        s += ", ";
    }
    s += "), {";
    for (size_t i = 0; i < defining_points.size(); ++i) {
      s += defining_points[i];
      if (i != defining_points.size() - 1)
        s += ", ";
    }
    s += "}";
    return s;
  }
};

struct Solution {
 public:
  std::string id;
  hv_t hypervolume = 0;
  hv_t real_hypervolume = 0;
  std::vector<int64_t> coordinates;
  std::vector<RefPoint> dominated_bounds;

  Solution(std::string id = "z0", hv_t hypervolume = 0, hv_t real_hypervolume = 0, std::vector<int64_t> coordinates = {}, std::vector<RefPoint> dominated_bounds = {})
      : id(std::move(id)), hypervolume(hypervolume), real_hypervolume(real_hypervolume), coordinates(std::move(coordinates)), dominated_bounds(std::move(dominated_bounds)) {}

  bool operator==(const Solution& other) const {
    return coordinates == other.coordinates;
  }

  std::string to_string() const {
    std::string s = id + " (";
    for (size_t i = 0; i < coordinates.size(); ++i) {
      s += std::to_string(coordinates[i]);
      if (i != coordinates.size() - 1)
        s += ", ";
    }
    s += ")";
    if (hypervolume != -1) {
      s += ", " + moco::ghps::to_string(hypervolume);
    }
    if (!dominated_bounds.empty()) {
      s += ", " + std::to_string(dominated_bounds.size());
    }
    s += "\n";
    return s;
  }

  static std::vector<Solution> convertSolutions(const std::vector<std::vector<int64_t>>& coordinates) {
    std::vector<Solution> solutions;
    solutions.reserve(coordinates.size());
    for (size_t i = 0; i < coordinates.size(); ++i) {
      solutions.emplace_back("z" + std::to_string(i), 0, 0, coordinates[i], std::vector<RefPoint>());
    }
    return solutions;
  }
};

} // namespace moco::ghps

#endif  // GHPS_CORE_TYPES_HPP
