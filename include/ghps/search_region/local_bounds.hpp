#ifndef GHPS_SEARCH_REGION_LOCAL_BOUNDS_HPP
#define GHPS_SEARCH_REGION_LOCAL_BOUNDS_HPP

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <ghps/core/types.hpp>
#include <ghps/core/timer.hpp>

namespace moco::ghps {


// Local lower bounds for maximization; search-region representation background:
// Klamroth, Lacour, Vanderpooten (2015), https://arxiv.org/abs/1502.06111
class LLB {
 public:
  // Constructor
  LLB() = default;

  static void update_bounds(std::vector<RefPoint> &bounds, Solution &solution, bool validate = false) {
    if (validate) {
      auto local_bounds_naive = bounds;
      auto local_bounds_improved = bounds;
      update_bounds_naive(local_bounds_naive, solution);
      update_bounds_article(local_bounds_improved, solution);
      if (local_bounds_naive.size() != local_bounds_improved.size()) {
        std::cout << "Error: The number of bounds is different in the naive and improved implementations" << std::endl;
        std::cout << "Naive: " << local_bounds_naive.size() << std::endl;
        for (const auto &bound : local_bounds_naive) {
          std::cout << bound.to_string() << std::endl;
        }
        std::cout << "Improved: " << local_bounds_improved.size() << std::endl;
        for (const auto &bound : local_bounds_improved) {
          std::cout << bound.to_string() << std::endl;
        }
        throw std::logic_error("LLB validation failed: naive has " + std::to_string(local_bounds_naive.size()) + " bounds, improved has " + std::to_string(local_bounds_improved.size()));
      }
      std::sort(local_bounds_naive.begin(), local_bounds_naive.end(), sort_by_id);
      std::sort(local_bounds_improved.begin(), local_bounds_improved.end(), sort_by_id);
      bool wrong_bounds = false;
      for (size_t i = 0; i < local_bounds_naive.size(); ++i) {
        if (local_bounds_naive[i].coordinates != local_bounds_improved[i].coordinates) {
          std::cout << "Error: The bounds are different in the naive and improved implementations" << std::endl;
          std::cout << "Naive: " << local_bounds_naive[i].to_string() << std::endl;
          std::cout << "Improved: " << local_bounds_improved[i].to_string() << std::endl;
          wrong_bounds = true;
        }
      }
      if (wrong_bounds) {
        throw std::logic_error("LLB validation failed: bounds differ between naive and improved implementations");
      }
    }
    update_bounds_article(bounds, solution);
  }

  // Minimization counterpart of the local-bound update.
  static void update_bounds_minimization(std::vector<RefPoint> &bounds, Solution &solution, bool validate = false) {
    if (validate) {
      auto local_bounds_naive = bounds;
      auto local_bounds_improved = bounds;
      update_bounds_naive_min(local_bounds_naive, solution);
      update_bounds_article_min(local_bounds_improved, solution);
      if (local_bounds_naive.size() != local_bounds_improved.size()) {
        std::cout << "Error: The number of bounds is different in the naive and improved implementations" << std::endl;
        std::cout << "Naive: " << local_bounds_naive.size() << std::endl;
        for (const auto &bound : local_bounds_naive) {
          std::cout << bound.to_string() << std::endl;
        }
        std::cout << "Improved: " << local_bounds_improved.size() << std::endl;
        for (const auto &bound : local_bounds_improved) {
          std::cout << bound.to_string() << std::endl;
        }
        throw std::logic_error("LLB min validation failed: naive has " + std::to_string(local_bounds_naive.size()) + " bounds, improved has " + std::to_string(local_bounds_improved.size()));
      }
      std::sort(local_bounds_naive.begin(), local_bounds_naive.end(), sort_by_id);
      std::sort(local_bounds_improved.begin(), local_bounds_improved.end(), sort_by_id);
      bool wrong_bounds = false;
      for (size_t i = 0; i < local_bounds_naive.size(); ++i) {
        if (local_bounds_naive[i].coordinates != local_bounds_improved[i].coordinates) {
          std::cout << "Error: The bounds are different in the naive and improved implementations" << std::endl;
          std::cout << "Naive: " << local_bounds_naive[i].to_string() << std::endl;
          std::cout << "Improved: " << local_bounds_improved[i].to_string() << std::endl;
          wrong_bounds = true;
        }
      }
      if (wrong_bounds) {
        throw std::logic_error("LLB min validation failed: bounds differ between naive and improved implementations");
      }
    }
    update_bounds_article_min(bounds, solution);
  }

  // Get dominated bounds
  static std::vector<RefPoint> get_dominated_bounds(const std::vector<RefPoint> &bounds, const Solution &sol) {
    std::vector<RefPoint> dominated_bounds;
    for (const auto &bound : bounds) {
      if (strongly_dominates(sol.coordinates, bound.coordinates)) {
        dominated_bounds.push_back(bound);
      }
    }
    return dominated_bounds;
  }

  // Get dominated bounds
  static std::vector<int32_t> get_dominated_bounds_idx(const std::vector<RefPoint> &bounds, const Solution &sol) {
    std::vector<int32_t> dominated_bounds;
    for (int i = 0; i < (int)bounds.size(); i++) {
      if (strongly_dominates(sol.coordinates, bounds[i].coordinates)) {
        dominated_bounds.push_back(i);
      }
    }
    return dominated_bounds;
  }

  template <bool IsMaximization = true>
  static bool weakly_dominates_impl(const std::vector<int64_t> &v1, const std::vector<int64_t> &v2) {
    for (size_t i = 0; i < v1.size(); ++i) {
      if constexpr (IsMaximization) {
        if (v1[i] < v2[i]) return false;
      } else {
        if (v1[i] > v2[i]) return false;
      }
    }
    return true;
  }

  template <bool IsMaximization = true>
  static bool strongly_dominates_impl(const std::vector<int64_t> &v1, const std::vector<int64_t> &v2) {
    for (size_t i = 0; i < v1.size(); ++i) {
      if constexpr (IsMaximization) {
        if (v1[i] <= v2[i]) return false;
      } else {
        if (v1[i] >= v2[i]) return false;
      }
    }
    return true;
  }

  template <bool IsMaximization = true>
  static std::vector<RefPoint> create_bounds(const RefPoint &ref, const Solution &sol) {
    std::vector<RefPoint> bounds;
    if (!strongly_dominates_impl<IsMaximization>(sol.coordinates, ref.coordinates)) {
      return bounds;
    }
    assert(ref.coordinates.size() == sol.coordinates.size());
    int32_t dimensions = ref.coordinates.size();
    for (int32_t i = 0; i < dimensions; i++) {
      std::string bound_id = ref.id + "." + std::to_string(i + 1);
      RefPoint new_bound = RefPoint(bound_id, ref.coordinates, ref.defining_points);
      new_bound.coordinates[i] = sol.coordinates[i];
      new_bound.defining_points[i] = sol.id;
      bounds.push_back(new_bound);
    }
    return bounds;
  }

 public:
  static bool weakly_dominates(const std::vector<int64_t> &v1, const std::vector<int64_t> &v2) {
    return weakly_dominates_impl<true>(v1, v2);
  }

 private:
  static bool strongly_dominates(const std::vector<int64_t> &v1, const std::vector<int64_t> &v2) {
    return strongly_dominates_impl<true>(v1, v2);
  }

  // Order bounds by identifier for validation.
  static bool sort_by_id(const RefPoint &a, const RefPoint &b) {
    assert(a.coordinates.size() == b.coordinates.size());
    return a.id < b.id;
  }

  template <bool IsMaximization = true>
  static void update_bounds_naive_impl(std::vector<RefPoint> &bounds, const Solution &sol) {
    std::vector<RefPoint> dominated_bounds;
    std::vector<int32_t> dominated_bounds_idx;
    for (int i = 0; i < (int)bounds.size(); i++) {
      if (strongly_dominates_impl<IsMaximization>(sol.coordinates, bounds[i].coordinates)) {
        dominated_bounds.push_back(bounds[i]);
        dominated_bounds_idx.push_back(i);
      }
    }
    if (dominated_bounds.empty()) return;

    std::vector<RefPoint> new_bounds;
    for (auto &bound : dominated_bounds) {
      std::vector<RefPoint> aux_bounds = create_bounds<IsMaximization>(bound, sol);
      new_bounds.insert(new_bounds.end(), aux_bounds.begin(), aux_bounds.end());
    }

    for (auto it = dominated_bounds_idx.rbegin(); it != dominated_bounds_idx.rend(); it++) {
      std::swap(bounds[*it], bounds.back());
      bounds.pop_back();
    }

    std::vector<int32_t> redundant_bounds_idx;
    for (int i = 0; i < (int)new_bounds.size(); i++) {
      for (int j = i + 1; j < (int)new_bounds.size(); j++) {
        if constexpr (IsMaximization) {
          if (weakly_dominates_impl<true>(new_bounds[i].coordinates, new_bounds[j].coordinates)) {
            redundant_bounds_idx.push_back(i);
          } else if (weakly_dominates_impl<true>(new_bounds[j].coordinates, new_bounds[i].coordinates)) {
            redundant_bounds_idx.push_back(j);
          }
        } else {
          if (weakly_dominates_impl<false>(new_bounds[i].coordinates, new_bounds[j].coordinates)) {
            redundant_bounds_idx.push_back(j);
          } else if (weakly_dominates_impl<false>(new_bounds[j].coordinates, new_bounds[i].coordinates)) {
            redundant_bounds_idx.push_back(i);
          }
        }
      }
    }

    std::sort(redundant_bounds_idx.begin(), redundant_bounds_idx.end());
    redundant_bounds_idx.erase(std::unique(redundant_bounds_idx.begin(), redundant_bounds_idx.end()), redundant_bounds_idx.end());
    for (auto it = redundant_bounds_idx.rbegin(); it != redundant_bounds_idx.rend(); ++it) {
      std::swap(new_bounds[*it], new_bounds.back());
      new_bounds.pop_back();
    }

    redundant_bounds_idx.clear();
    for (int j = 0; j < (int)bounds.size(); j++) {
      for (int i = 0; i < (int)new_bounds.size(); i++) {
        if constexpr (IsMaximization) {
          if (weakly_dominates_impl<true>(new_bounds[i].coordinates, bounds[j].coordinates)) {
            redundant_bounds_idx.push_back(i);
          }
        } else {
          if (weakly_dominates_impl<false>(bounds[j].coordinates, new_bounds[i].coordinates)) {
            redundant_bounds_idx.push_back(i);
          }
        }
      }
    }

    std::sort(redundant_bounds_idx.begin(), redundant_bounds_idx.end());
    redundant_bounds_idx.erase(std::unique(redundant_bounds_idx.begin(), redundant_bounds_idx.end()), redundant_bounds_idx.end());
    for (auto it = redundant_bounds_idx.rbegin(); it != redundant_bounds_idx.rend(); ++it) {
      std::swap(new_bounds[*it], new_bounds.back());
      new_bounds.pop_back();
    }

    bounds.insert(bounds.end(), new_bounds.begin(), new_bounds.end());
  }

  template <bool IsMaximization = true>
  static void update_bounds_article_impl(std::vector<RefPoint> &bounds, Solution &sol, const Timer& deadline) {
    deadline.check();
    int32_t dimensions = sol.coordinates.size();
    std::vector<RefPoint> A;
    std::vector<int64_t> A_idx;
    for (int i = 0; i < (int)bounds.size(); i++) {
      deadline.check();
      if (strongly_dominates_impl<IsMaximization>(sol.coordinates, bounds[i].coordinates)) {
        A.push_back(bounds[i]);
        A_idx.push_back(i);
      }
    }
    if (A.empty()) return;

    std::vector<std::vector<RefPoint>> B(dimensions, std::vector<RefPoint>());
    for (int j = 0; j < dimensions; j++) {
      deadline.check();
      for (const auto &bound : bounds) {
        deadline.check();
        if (sol.coordinates[j] == bound.coordinates[j]) {
          bool is_dominated = false;
          for (int k = 0; k < dimensions; k++) {
            deadline.check();
            if (k != j) {
              if constexpr (IsMaximization) {
                if (sol.coordinates[k] < bound.coordinates[k]) {
                  is_dominated = true;
                  break;
                }
              } else {
                if (sol.coordinates[k] > bound.coordinates[k]) {
                  is_dominated = true;
                  break;
                }
              }
            }
          }
          if (!is_dominated) {
            B[j].push_back(bound);
          }
        }
      }
    }

    std::vector<std::vector<RefPoint>> P(dimensions, std::vector<RefPoint>());
    for (const auto &bound : A) {
      deadline.check();
      for (int i = 0; i < dimensions; i++) {
        deadline.check();
        std::string bound_id = bound.id + "." + std::to_string(i + 1);
        RefPoint new_bound = RefPoint(bound_id, bound.coordinates, bound.defining_points);
        new_bound.coordinates[i] = sol.coordinates[i];
        new_bound.defining_points[i] = sol.id;
        P[i].push_back(new_bound);
      }
    }

    for (int j = 0; j < dimensions; j++) {
      deadline.check();
      std::vector<RefPoint> filtered_bounds;
      for (const auto &p : P[j]) {
        deadline.check();
        bool is_dominated = false;
        for (const auto &u : P[j]) {
          deadline.check();
          if (u.id != p.id) {
            if constexpr (IsMaximization) {
              if (weakly_dominates_impl<true>(p.coordinates, u.coordinates)) {
                is_dominated = true;
                break;
              }
            } else {
              if (weakly_dominates_impl<false>(u.coordinates, p.coordinates)) {
                is_dominated = true;
                break;
              }
            }
          }
        }
        if (is_dominated) continue;
        for (const auto &b : B[j]) {
          deadline.check();
          if constexpr (IsMaximization) {
            if (weakly_dominates_impl<true>(p.coordinates, b.coordinates)) {
              is_dominated = true;
              break;
            }
          } else {
            if (weakly_dominates_impl<false>(b.coordinates, p.coordinates)) {
              is_dominated = true;
              break;
            }
          }
        }
        if (!is_dominated) {
          filtered_bounds.push_back(p);
        }
      }
      P[j] = filtered_bounds;
    }

    for (auto it = A_idx.rbegin(); it != A_idx.rend(); ++it) {
      deadline.check();
      std::swap(bounds[*it], bounds.back());
      bounds.pop_back();
    }

    for (int i = 0; i < dimensions; i++) {
      deadline.check();
      bounds.insert(bounds.end(), P[i].begin(), P[i].end());
    }
  }

 public:
  static void update_bounds_naive(std::vector<RefPoint> &bounds, const Solution &sol) {
    update_bounds_naive_impl<true>(bounds, sol);
  }
  static void update_bounds_article(std::vector<RefPoint> &bounds, Solution &sol, const Timer& deadline = Timer()) {
    update_bounds_article_impl<true>(bounds, sol, deadline);
  }
  static void update_bounds_naive_min(std::vector<RefPoint> &bounds, const Solution &sol) {
    update_bounds_naive_impl<false>(bounds, sol);
  }
  static void update_bounds_article_min(std::vector<RefPoint> &bounds, Solution &sol, const Timer& deadline = Timer()) {
    update_bounds_article_impl<false>(bounds, sol, deadline);
  }
};

} // namespace moco::ghps

#endif  // GHPS_SEARCH_REGION_LOCAL_BOUNDS_HPP
