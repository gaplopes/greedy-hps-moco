#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <ghps/core/hypervolume.hpp>
#include <ghps/core/nondominated_data.hpp>
#include <ghps/core/types.hpp>
#include <ghps/mokp/instance.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace moco::ghps {

class Statistics {
 private:
  const std::shared_ptr<const detail::NondominatedData> front_data;

 public:
  const std::vector<std::vector<int64_t>>& nondominated_set;
  std::vector<Solution> solutions;
  std::vector<std::tuple<int32_t, int64_t, double>> statistics;
  double elapsed_time = 0;
  int64_t iterations = 0;
  int64_t n_subproblems = 0;
  int64_t n_cache_hits = 0;
  double subgraphs_time = 0;
  double subproblems_time = 0;
  const bool is_maximization;
  SolverStatus status = SolverStatus::Optimal;
  std::string status_message;

  Statistics(const MOKP& problem)
      : Statistics(problem, {}, {}, 0, 0, 0, 0, 0, 0) {}

  Statistics(const MOKP& problem,
             std::vector<Solution> solutions,
             std::vector<std::tuple<int32_t, int64_t, double>> statistics,
             double elapsed_time,
             int64_t iterations,
             int64_t n_subproblems,
             int64_t n_cache_hits,
             double subgraphs_time,
             double subproblems_time)
      : front_data(problem.reporting_data()), nondominated_set(front_data->points),
        solutions(std::move(solutions)), statistics(std::move(statistics)),
        elapsed_time(elapsed_time), iterations(iterations), n_subproblems(n_subproblems),
        n_cache_hits(n_cache_hits), subgraphs_time(subgraphs_time), subproblems_time(subproblems_time),
        is_maximization(front_data->is_maximization), reference_point(front_data->reference), nd_set_values(front_data->membership) {}

  std::string to_string(bool detailed = false) const {
    const auto selected = coordinates(solutions);
    const int n_matching = std::count_if(selected.begin(), selected.end(),
        [this](const auto& point) { return nd_set_values.contains(point); });
    const hv_t solutions_hv = calculate_hv(reference_point, selected);
    const hv_t nondominated_hv = front_data->hypervolume;
    const double ratio_hv = nondominated_hv > 0 ? double(solutions_hv) / double(nondominated_hv) : 0.0;
    const auto [nadir_nset_hv, nadir_set_hv, ratio_nadir_hv] = get_nadir_hvs();
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << "Status: " << status_name(status) << '\n';
    if (!status_message.empty()) out << status_message << '\n';
    int64_t total_recursions = 0;
    if (!statistics.empty()) {
      if (detailed) out << "Statistics:\n";
      for (size_t i = 0; i < statistics.size(); ++i) {
        const auto [k, recursions, time] = statistics[i];
        if (i > 0) out << ' ';
        out << '(' << k << ',' << recursions << ',' << time << ')';
        total_recursions += recursions;
      }
      out << '\n';
    }
    if (detailed) {
      out << "Iterations: " << iterations << '\n'
          << "Total recursions: " << total_recursions << '\n'
          << "Subproblems: " << n_subproblems << '\n'
          << "Cache hits: " << n_cache_hits << '\n'
          << "Subgraphs time: " << subgraphs_time << "s\n"
          << "Subproblems time: " << subproblems_time << "s\n";
    } else {
      out << iterations << ' ' << n_subproblems << ' ' << n_cache_hits << ' '
          << subgraphs_time << ' ' << subproblems_time << '\n';
    }
    if (!solutions.empty()) {
      if (detailed) out << "Solutions:\n";
      for (size_t i = 0; i < solutions.size(); ++i) {
        if (i > 0) out << ' ';
        out << '(';
        for (size_t j = 0; j < solutions[i].coordinates.size(); ++j) {
          if (j > 0) out << ',';
          out << solutions[i].coordinates[j];
        }
        out << ")[" << solutions[i].hypervolume << ']';
      }
      out << '\n';
    }
    if (detailed) {
      out << "Elapsed time: " << elapsed_time << "s\n"
          << "|Nset|: " << nondominated_set.size() << '\n'
          << "|Solutions|: " << solutions.size() << '\n'
          << "|Matching|: " << n_matching << '\n'
          << "HV(Nset): " << nondominated_hv << '\n'
          << "HV(Solutions): " << solutions_hv << '\n'
          << "HV ratio: " << ratio_hv << '\n'
          << "HV(Nadir Nset): " << nadir_nset_hv << '\n'
          << "HV(Nadir Solutions): " << nadir_set_hv << '\n'
          << "HV ratio (Nadir): " << ratio_nadir_hv;
    } else {
      out << nondominated_hv << ' ' << solutions_hv << ' ' << ratio_hv << '\n'
          << nadir_nset_hv << ' ' << nadir_set_hv << ' ' << ratio_nadir_hv << '\n'
          << nondominated_set.size() << ' ' << solutions.size() << ' ' << n_matching << ' ' << elapsed_time;
    }
    return out.str();
  }

  void save_sequence(const std::filesystem::path& path) const {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Could not open file " + path.string() + " for writing");
    for (const auto& sol : solutions) {
      for (int64_t value : sol.coordinates) out << value << ' ';
      out << sol.hypervolume << '\n';
    }
    const auto [nset_hv, solutions_hv, ratio] = get_nadir_hvs();
    out << nset_hv << ' ' << solutions_hv << ' ' << ratio << '\n';
    out.close();
    if (!out) throw std::runtime_error("Could not write file " + path.string());
  }

  void compare_solution_sets(const std::vector<Solution>& set1, const std::vector<Solution>& set2) const {
    if (set1.empty() || set2.empty()) return;
    const auto points1 = coordinates(set1), points2 = coordinates(set2);
    const auto nadir = compute_nadir(points2);
    const hv_t hv1 = calculate_hv(reference_point, points1), hv2 = calculate_hv(reference_point, points2);
    const hv_t nadir_hv1 = calculate_hv(nadir, points1), nadir_hv2 = calculate_hv(nadir, points2);
    std::cout << hv1 << ' ' << hv2 << ' ' << (hv2 > 0 ? double(hv1) / double(hv2) : 0.0) << '\n'
              << nadir_hv1 << ' ' << nadir_hv2 << ' ' << (nadir_hv2 > 0 ? double(nadir_hv1) / double(nadir_hv2) : 0.0) << '\n'
              << set1.size() << ' ' << set2.size() << ' ' << matching(points1, points2) << '\n';
  }

  std::tuple<hv_t, hv_t, double> get_nadir_hvs() const {
    if (nondominated_set.empty()) return {0, 0, 0};
    const hv_t nset_hv = front_data->nadir_hypervolume;
    const hv_t solutions_hv = calculate_hv(front_data->nadir, coordinates(solutions));
    return {nset_hv, solutions_hv, nset_hv > 0 ? double(solutions_hv) / double(nset_hv) : 0.0};
  }

 private:
  const std::vector<int64_t>& reference_point;
  const std::set<std::vector<int64_t>>& nd_set_values;

  static std::vector<std::vector<int64_t>> coordinates(const std::vector<Solution>& solutions) {
    std::vector<std::vector<int64_t>> result;
    result.reserve(solutions.size());
    for (const auto& sol : solutions) result.push_back(sol.coordinates);
    return result;
  }

  static int matching(const std::vector<std::vector<int64_t>>& reference,
                      const std::vector<std::vector<int64_t>>& target) {
    const std::set<std::vector<int64_t>> points(reference.begin(), reference.end());
    return std::count_if(target.begin(), target.end(), [&points](const auto& point) { return points.contains(point); });
  }

  std::vector<int64_t> compute_nadir(const std::vector<std::vector<int64_t>>& points) const {
    if (points.empty()) return reference_point;
    auto result = points.front();
    for (const auto& point : points) {
      for (size_t i = 0; i < point.size(); ++i)
        result[i] = is_maximization ? std::min(result[i], point[i]) : std::max(result[i], point[i]);
    }
    // Apply epsilon difference to ensure the nadir point is dominated by all solutions in the nondominated set
    for (auto& value : result) value += is_maximization ? -1 : 1;
    return result;
  }

  hv_t calculate_hv(const std::vector<int64_t>& reference,
                    const std::vector<std::vector<int64_t>>& points) const {
    std::vector<std::vector<hv_t>> wide_points;
    wide_points.reserve(points.size());
    for (const auto& point : points) wide_points.emplace_back(point.begin(), point.end());
    return HypervolumeIndicator<hv_t, std::vector<hv_t>>(
        std::vector<hv_t>(reference.begin(), reference.end()), is_maximization).set_hypervolume_unchecked(wide_points);
  }
};

} // namespace moco::ghps

#endif // STATISTICS_HPP
