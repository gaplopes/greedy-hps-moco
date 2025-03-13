#ifndef STRUCTS_HPP_
#define STRUCTS_HPP_

#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

struct RefPoint {
 public:
  std::string id;
  std::vector<int64_t> coordinates;
  std::vector<std::string> defining_points;

  RefPoint() = default;

  RefPoint(const std::string& id, const std::vector<int64_t>& coordinates, const std::vector<std::string>& defining_points = {})
      : id(id), coordinates(coordinates), defining_points(defining_points) {}

  bool operator==(const RefPoint& other) const {
    return coordinates == other.coordinates;
  }

  static RefPoint origin(const int32_t &M, std::vector<int64_t> coordinates = {}) {
    if (coordinates.empty()) {
      coordinates = std::vector<int64_t>(M, 0);
    }
    std::string id = "u0";
    std::vector<std::string> defining_points;
    for (int i = 0; i < M; i++) {
      defining_points.push_back("d" + std::to_string(i + 1));
    }
    return RefPoint(id, coordinates, defining_points);
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
  int64_t hypervolume;
  int64_t real_hypervolume;
  std::vector<int64_t> coordinates;
  std::vector<RefPoint> dominated_bounds;

  Solution(const std::string& id = "z0", int64_t hypervolume = 0, int64_t real_hypervolume = 0, const std::vector<int64_t>& coordinates = {}, const std::vector<RefPoint>& dominated_bounds = {})
      : id(id), hypervolume(hypervolume), real_hypervolume(real_hypervolume), coordinates(coordinates), dominated_bounds(dominated_bounds) {}

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
      s += ", " + std::to_string(hypervolume);
    }
    if (!dominated_bounds.empty()) {
      s += ", " + std::to_string(dominated_bounds.size());
    }
    s += "\n";
    return s;
  }

  static std::vector<Solution> convertSolutions(const std::vector<std::vector<int64_t>>& coordinates) {
    std::vector<Solution> solutions;
    for (size_t i = 0; i < coordinates.size(); ++i) {
      solutions.emplace_back("z" + std::to_string(i), 0, 0, coordinates[i], std::vector<RefPoint>());
    }
    return solutions;
  }
};

#endif  // STRUCTS_HPP_
