#ifndef GHPS_EXPERIMENT_HELPER_HPP
#define GHPS_EXPERIMENT_HELPER_HPP

#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <ghps/core/hypervolume.hpp>
#include <ghps/core/types.hpp>

#include "config.hpp"

namespace moco::ghps {

inline const std::map<Solver, const char*> solver_csv_filenames = {
    {Solver::Milp, "results_milp.csv"}, {Solver::Miqp, "results_miqp.csv"},
    {Solver::Log, "results_log.csv"}, {Solver::Greedy, "results_greedy.csv"},
    {Solver::BranchBound, "results_bb.csv"}};

inline const std::string csv_header =
    "m,n,seed,k,refs,nset_hv,diff_v,sim_sol,sim_hv,sol,"
    "hv,hv_real,time,nodes,status,backend_status,message\n";

// Auxiliary functions

template <typename T>
inline std::string v_to_str(const std::vector<T> &coords) {
  std::string str = "[";
  str.reserve(coords.size() * 24 + 2);
  for (size_t i = 0; i < coords.size(); ++i) {
    str += std::to_string(coords[i]);
    if (i < coords.size() - 1) {
      str += "; ";
    }
  }
  str += "]";
  return str;
}

inline std::string v_to_str(const RefPoint &vec_coord) {
  return v_to_str(vec_coord.coordinates);
}

inline std::string v_to_str(const std::vector<RefPoint> &vec_coords) {
  std::string str = "[";
  for (size_t i = 0; i < vec_coords.size(); ++i) {
    str += v_to_str(vec_coords[i]);
    if (i < vec_coords.size() - 1) {
      str += "; ";
    }
  }
  str += "]";
  return str;
}

inline std::string diff_v(const std::vector<int64_t> &r,
                          const std::vector<int64_t> &s) {
  std::vector<int64_t> result(r.size());
  for (size_t i = 0; i < r.size(); ++i) {
    result[i] = std::abs(s[i] - r[i]);
  }
  return v_to_str(result);
}

inline void
write_line_to_csv(std::ofstream &csv_file, int m, int n,
                     int seed, int k, const std::string &refs, hv_t nset_hv,
                     const std::string &diff_v_str, const std::string &sim_sol,
                     hv_t sim_hv, const std::string &sol, hv_t hv,
                     hv_t hv_real, double time, int64_t nodes = 0, SolverStatus status = SolverStatus::Optimal,
                     int backend_status = 0, const std::string& message = {}) {
  csv_file << m << "," << n << "," << seed << "," << k << "," << refs << ","
           << nset_hv << "," << diff_v_str << "," << sim_sol << "," << sim_hv
           << "," << sol << "," << hv << "," << hv_real << "," << time << ","
           << nodes << "," << status_name(status) << "," << backend_status << ",\"";
  for (char c : message) {
    if (c == '"') csv_file << '"';
    csv_file << (c == '\n' || c == '\r' ? ' ' : c);
  }
  csv_file << "\"\n";
  csv_file.flush();
  if (!csv_file) throw std::runtime_error("Could not write CSV results");
}

inline void init_csv_file(std::ofstream &csv_file, const std::string &csv_filename,
                   const std::string &csv_header) {
  // Open in append mode; write header only if file is missing or empty
  bool write_header = true;
  {
    std::ifstream check(csv_filename);
    write_header = !check.good() || check.peek() == std::ifstream::traits_type::eof();
    if (!write_header) {
      std::string header;
      std::getline(check, header);
      if (header + "\n" != csv_header)
        throw std::runtime_error("CSV schema differs; use a fresh results directory: " + csv_filename);
    }
  }

  csv_file.open(csv_filename, std::ios::app);
  if (!csv_file.is_open()) {
    throw std::runtime_error("Could not open CSV file for writing: " + csv_filename);
  }
  if (write_header) {
    csv_file << csv_header;
    csv_file.flush();
    if (!csv_file) throw std::runtime_error("Could not write CSV header: " + csv_filename);
  }
}

inline void compute_hypervolume_from_strings(const std::string &refs_str,
                                             const std::string &sol_str) {
  auto parse_matrix = [](const std::string &str) {
    std::vector<std::vector<int64_t>> points;
    std::vector<int64_t> coords;
    std::string num_str;
    for (char c : str) {
      if (std::isdigit(static_cast<unsigned char>(c)) || c == '-') {
        num_str += c;
      } else {
        if (!num_str.empty()) {
          coords.push_back(std::stoll(num_str));
          num_str.clear();
        }
        if (c == ']' && !coords.empty()) {
          points.push_back(coords);
          coords.clear();
        }
      }
    }
    return points;
  };

  auto parse_vector = [](const std::string &str) {
    std::vector<int64_t> coords;
    std::string num_str;
    for (char c : str) {
      if (std::isdigit(static_cast<unsigned char>(c)) || c == '-') {
        num_str += c;
      } else {
        if (!num_str.empty()) {
          coords.push_back(std::stoll(num_str));
          num_str.clear();
        }
      }
    }
    return coords;
  };

  std::vector<std::vector<int64_t>> ref_points = parse_matrix(refs_str);
  std::vector<int64_t> sol = parse_vector(sol_str);

  hv_t hv = scalarized_hypervolume(ref_points, sol);
  std::cout << "Hypervolume: " << hv << std::endl;
}

} // namespace moco::ghps

#endif // GHPS_EXPERIMENT_HELPER_HPP
