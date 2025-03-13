#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <hvi.hpp>
#include <mokp.hpp>
#include <string>
#include <structs.hpp>
#include <tuple>
#include <vector>

class Statistics {
 public:
  MOKP problem;
  std::vector<Solution> nondominated_set;
  std::vector<Solution> solutions;
  std::vector<std::tuple<int32_t, int64_t, double>> statistics;
  double elapsed_time;
  int64_t iterations = 0;
  int64_t n_subproblems = 0;
  int64_t n_cache_hits = 0;
  double subgraphs_time = 0;
  double subproblems_time = 0;
  bool is_maximization = true;

  Statistics(MOKP problem) : problem(problem) {}

  Statistics(MOKP problem,
             std::vector<Solution> solutions,
             std::vector<std::tuple<int32_t, int64_t, double>> statistics,
             double elapsed_time,
             int64_t iterations,
             int64_t n_subproblems,
             int64_t n_cache_hits,
             double subgraphs_time,
             double subproblems_time,
             bool is_maximization = true) : problem(problem),
                                            solutions(solutions),
                                            statistics(statistics),
                                            elapsed_time(elapsed_time),
                                            iterations(iterations),
                                            n_subproblems(n_subproblems),
                                            n_cache_hits(n_cache_hits),
                                            subgraphs_time(subgraphs_time),
                                            subproblems_time(subproblems_time),
                                            is_maximization(is_maximization) {
    solutions_set_hv = calculate_hv(problem.reference_point, this->solutions);
    n_solutions = this->solutions.size();
    nondominated_set = problem.nondominated_set;
    if (!nondominated_set.empty()) {
      n_nondominated_set = nondominated_set.size();
      // Nset vs solutions
      nondominated_set_hv = calculate_hv(problem.reference_point, nondominated_set);
      ratio_hv = static_cast<double>(solutions_set_hv) / nondominated_set_hv;
      // Nset vs solutions (using nadir)
      std::vector<int64_t> nadir_nset = compute_nadir_nset(nondominated_set);
      nadir_nset_hv = calculate_hv(nadir_nset, nondominated_set);
      nadir_set_hv = calculate_hv(nadir_nset, this->solutions);
      ratio_nadir_hv = static_cast<double>(nadir_set_hv) / nadir_nset_hv;
      // Matching
      n_matching = matching(nondominated_set, this->solutions);
    }
  }

  std::string to_string(bool detailed = false) const {
    std::string stats_str;

    if (!statistics.empty()) {
      if (detailed) {
        stats_str += "Statistics:\n";
      }
      for (const auto& stat : statistics) {
        stats_str += "(" + std::to_string(std::get<0>(stat)) + "," + std::to_string(std::get<1>(stat)) + "," + std::to_string(std::get<2>(stat)) + ") ";
      }
      stats_str.pop_back();  // Remove the trailing space
      stats_str += "\n";
    }

    if (detailed) {
      stats_str += "Iterations: " + std::to_string(iterations) + "\n";
      stats_str += "Subproblems: " + std::to_string(n_subproblems) + "\n";
      stats_str += "Cache hits: " + std::to_string(n_cache_hits) + "\n";
      stats_str += "Subgraphs time: " + std::to_string(subgraphs_time) + "s\n";
      stats_str += "Subproblems time: " + std::to_string(subproblems_time) + "s\n";
    } else {
      stats_str += std::to_string(iterations) + " ";
      stats_str += std::to_string(n_subproblems) + " ";
      stats_str += std::to_string(n_cache_hits) + " ";
      stats_str += std::to_string(subgraphs_time) + " ";
      stats_str += std::to_string(subproblems_time) + "\n";
    }

    if (!solutions.empty()) {
      if (detailed) {
        stats_str += "Solutions:\n";
      }
      for (const auto& sol : solutions) {
        stats_str += "(";
        for (std::size_t j = 0; j < sol.coordinates.size(); ++j) {
          stats_str += std::to_string(sol.coordinates[j]);
          if (j < sol.coordinates.size() - 1) {
            stats_str += ",";
          }
        }
        stats_str += ") ";
      }
      stats_str.pop_back();  // Remove the trailing space
      stats_str += "\n";
    }

    if (detailed) {
      stats_str += "Elapsed time: " + std::to_string(elapsed_time) + "s\n";
      stats_str += "|Nset|: " + std::to_string(n_nondominated_set) + "\n";
      stats_str += "|Solutions|: " + std::to_string(n_solutions) + "\n";
      stats_str += "|Matching|: " + std::to_string(n_matching) + "\n";
      stats_str += "HV(Nset): " + std::to_string(nondominated_set_hv) + "\n";
      stats_str += "HV(Solutions): " + std::to_string(solutions_set_hv) + "\n";
      stats_str += "HV ratio: " + std::to_string(ratio_hv) + "\n";
      stats_str += "HV(Nadir Nset): " + std::to_string(nadir_nset_hv) + "\n";
      stats_str += "HV(Nadir Solutions): " + std::to_string(nadir_set_hv) + "\n";
      stats_str += "HV ratio (Nadir): " + std::to_string(ratio_nadir_hv);
    } else {
      stats_str += std::to_string(nondominated_set_hv) + " " + std::to_string(solutions_set_hv) + " " + std::to_string(ratio_hv) + "\n";
      stats_str += std::to_string(nadir_nset_hv) + " " + std::to_string(nadir_set_hv) + " " + std::to_string(ratio_nadir_hv) + "\n";
      stats_str += std::to_string(n_nondominated_set) + " " + std::to_string(n_solutions) + " ";
      stats_str += std::to_string(n_matching) + " " + std::to_string(elapsed_time);
    }
    return stats_str;
  }

  void compare_solution_sets(const std::vector<Solution>& set1, const std::vector<Solution>& set2) const {
    // Check if the sets are empty
    if (set1.empty() || set2.empty()) {
      return;
    }
    // Compare the hypervolume of the two sets
    int64_t hv_set1 = calculate_hv(problem.reference_point, set1);
    int64_t hv_set2 = calculate_hv(problem.reference_point, set2);
    double ratio = static_cast<double>(hv_set1) / hv_set2;
    // Compare the hypervolume of the two sets (using nadir)
    // Note: The nadir set is computed using the set2
    std::vector<int64_t> nadir_point = compute_nadir_nset(set2);
    int64_t nadir_hv_set1 = calculate_hv(nadir_point, set1);
    int64_t nadir_hv_set2 = calculate_hv(nadir_point, set2);
    double ratio_nadir = static_cast<double>(nadir_hv_set1) / nadir_hv_set2;
    // Compare the number of solutions
    int n_set1 = set1.size();
    int n_set2 = set2.size();
    // Compare the number of matching solutions
    int n_matching = matching(set1, set2);
    // Print the results
    std::cout << hv_set1 << " " << hv_set2 << " " << ratio << std::endl;
    std::cout << nadir_hv_set1 << " " << nadir_hv_set2 << " " << ratio_nadir << std::endl;
    std::cout << n_set1 << " " << n_set2 << " " << n_matching << std::endl;
  }

 private:
  int n_solutions = 0;
  int64_t solutions_set_hv;

  int n_nondominated_set = 0;
  int n_matching = 0;

  int64_t nondominated_set_hv = 0;
  double ratio_hv = 0;

  int64_t nadir_nset_hv = 0;
  int64_t nadir_set_hv = 0;
  double ratio_nadir_hv = 0;

  int matching(const std::vector<Solution>& nondominated_set, const std::vector<Solution>& solutions) const {
    if (nondominated_set.empty() || solutions.empty()) {
      return 0;
    }
    int matching = 0;
    for (const auto& sol : solutions) {
      if (std::find(nondominated_set.begin(), nondominated_set.end(), sol) != nondominated_set.end()) {
        matching++;
      }
    }
    return matching;
  }

  std::vector<int64_t> compute_nadir_nset(const std::vector<Solution>& nondominated_set) const {
    if (nondominated_set.empty()) {
      return problem.reference_point;
    }
    std::vector<int64_t> nadir(problem.M, 0);
    if (is_maximization) {
      for (auto& val : nadir) {
        val = std::numeric_limits<int64_t>::max();
      }
    } else {
      for (auto& val : nadir) {
        val = std::numeric_limits<int64_t>::min();
      }
    }
    for (const auto& sol : nondominated_set) {
      for (std::size_t i = 0; i < sol.coordinates.size(); ++i) {
        if (is_maximization) {
          nadir[i] = std::min(nadir[i], sol.coordinates[i]);
        } else {
          nadir[i] = std::max(nadir[i], sol.coordinates[i]);
        }
      }
    }
    return nadir;
  }

  int64_t calculate_hv(const std::vector<int64_t>& ref_point, const std::vector<Solution>& solutions) const {
    if (solutions.empty()) {
      return 0;
    }
    std::vector<std::vector<int64_t>> solutions_coordinates;
    for (const auto& sol : solutions) {
      solutions_coordinates.push_back(sol.coordinates);
    }
    return HypervolumeIndicator<int64_t, std::vector<int64_t>>(ref_point, is_maximization).set_hypervolume(solutions_coordinates);
  }
};

#endif  // STATISTICS_HPP
