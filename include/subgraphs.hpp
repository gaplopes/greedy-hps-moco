#ifndef SUBGRAPHS_HPP
#define SUBGRAPHS_HPP

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <rwd.hpp>
#include <set>
#include <sstream>
#include <string>
#include <utils.hpp>
#include <vector>

class Node {
 public:
  std::string id;
  std::vector<std::string> defining_points;
  std::vector<Node> neighbors;

  Node(std::string id, std::vector<std::string> defining_points)
      : id(std::move(id)), defining_points(std::move(defining_points)), neighbors() {}

  std::string to_string() const {
    std::string str = id + " {";
    str += defining_points[0];
    for (size_t i = 1; i < defining_points.size(); ++i) {
      str += ", " + defining_points[i];
    }
    str += "}, {";
    if (!neighbors.empty()) {
      str += neighbors[0].id;
      for (size_t i = 1; i < neighbors.size(); ++i) {
        str += ", " + neighbors[i].id;
      }
    }
    str += "}";
    return str;
  }
};

class Subgraphs {
 public:
  // Constructor
  Subgraphs() = default;

  // Get subgraphs from reference points
  static std::vector<std::vector<int32_t>> get_subgraphs(const std::vector<RefPoint>& local_bounds, int k, int time_limit = 120, bool validate = false) {
    if (validate && k > 1) {
      auto subsets_python = get_subgraph_enucon(local_bounds, k, 120);
      auto subsets_cpp = get_subgraphs_rwd(local_bounds, k, 120);
      if (subsets_python.size() != subsets_cpp.size()) {
        std::cout << "Error: The number of subsets is different in the Python and C++ implementations" << std::endl;
        std::cout << "Python: " << subsets_python.size() << std::endl;
        std::cout << "C++: " << subsets_cpp.size() << std::endl;
        std::cout << "Local bounds: " << local_bounds.size() << std::endl;
        exit(1);
      }
      // Sort the subsets
      for (auto& subset : subsets_python) {
        std::sort(subset.begin(), subset.end());
      }
      for (auto& subset : subsets_cpp) {
        std::sort(subset.begin(), subset.end());
      }
      std::sort(subsets_python.begin(), subsets_python.end());
      std::sort(subsets_cpp.begin(), subsets_cpp.end());
      if (subsets_python != subsets_cpp) {
        std::cout << "Error: The subsets are different in the Python and C++ implementations" << std::endl;
        exit(1);
      }
    }
    return get_subgraphs_rwd(local_bounds, k, time_limit);
  }

 private:
  // Check if two nodes are neighbors
  static bool are_neighbors(const Node& u, const Node& v) {
    assert(u.defining_points.size() == v.defining_points.size());
    size_t p = u.defining_points.size();
    std::set<size_t> indices;
    for (size_t i = 0; i < p; ++i) {
      indices.insert(i);
    }
    for (size_t i = 0; i < p; ++i) {
      if (u.defining_points[i] == v.defining_points[i]) {
        indices.erase(i);
      }
    }
    if (indices.size() != 2) {
      return false;
    }
    size_t j = *indices.begin();
    size_t k = *indices.rbegin();
    return (u.defining_points[j] == v.defining_points[k]) || (u.defining_points[k] == v.defining_points[j]);
  }

  // Get neighbors of a node
  static std::vector<Node> get_neighbors(const Node& u, const std::vector<Node>& nodes) {
    std::vector<Node> neighbors;
    for (const auto& v : nodes) {
      if (u.id != v.id && are_neighbors(u, v)) {
        neighbors.push_back(v);
      }
    }
    return neighbors;
  }

  // Compute the graph by finding neighbors for each node
  static void compute_graph(std::vector<Node>& nodes) {
    for (auto& u : nodes) {
      u.neighbors = get_neighbors(u, nodes);
    }
  }

  // Print nodes
  static void print_nodes(const std::vector<Node>& nodes) {
    for (const auto& n : nodes) {
      std::cout << n.to_string() << std::endl;
    }
  }

  // Print graph
  static void print_graph(const std::vector<Node>& nodes) {
    print_nodes(nodes);
  }

  // Get graph edges
  static std::vector<std::pair<int32_t, int32_t>> get_graph_edges(const std::vector<RefPoint>& ref_points) {
    std::vector<Node> nodes;
    for (const auto& ref : ref_points) {
      nodes.emplace_back(ref.id, ref.defining_points);
    }
    std::vector<std::pair<int32_t, int32_t>> edges;
    for (size_t i = 0; i < nodes.size(); ++i) {
      for (size_t j = i + 1; j < nodes.size(); ++j) {
        if (are_neighbors(nodes[i], nodes[j])) {
          edges.emplace_back(i, j);
        }
      }
    }
    return edges;
  }

  // Create graph file
  static void create_graph_file(const std::vector<RefPoint>& ref_points, const std::string& filename) {
    std::ofstream file(filename);
    auto edges = get_graph_edges(ref_points);
    for (const auto& e : edges) {
      file << e.first << " " << e.second << std::endl;
    }
  }

  // Create graph file from nodes
  static void create_graph_file(const std::vector<Node>& nodes, const std::string& filename) {
    std::ofstream file(filename);
    std::vector<std::pair<int32_t, int32_t>> edges;
    for (size_t i = 0; i < nodes.size(); ++i) {
      for (size_t j = i + 1; j < nodes.size(); ++j) {
        if (are_neighbors(nodes[i], nodes[j])) {
          edges.emplace_back(i, j);
        }
      }
    }
    for (const auto& e : edges) {
      file << e.first << " " << e.second << std::endl;
    }
  }

  // Get subgraphs from file
  static std::vector<std::vector<int32_t>> get_subgraph_enucon(const std::vector<RefPoint>& local_bounds, int k, int time_limit) {
    std::string file = "graph.in";
    create_graph_file(local_bounds, file);
    std::string command = "python3 ../others/enucon/enucon.py delay-old " +
                          std::to_string(time_limit) + " " + std::to_string(k) + " " + file +
                          " ../others/enucon/result/ > ../others/enucon/py_garbage.txt";
    int status = system(command.c_str());
    if (status == -1) {
      std::cerr << "Error while executing the command: " << command << std::endl;
      exit(1);
    }
    std::vector<std::vector<int32_t>> subgraphs;
    std::ifstream fin("py_subs.out");
    std::string line;
    while (std::getline(fin, line)) {
      std::vector<int32_t> subgraph;
      std::stringstream ss(line);
      std::string token;
      while (std::getline(ss, token, ' ')) {
        subgraph.push_back(std::stoi(token));
      }
      if (subgraph.size() == static_cast<size_t>(k)) {
        subgraphs.push_back(subgraph);
      }
    }
    return subgraphs;
  }

  // Get subgraphs from reference points
  static std::vector<std::vector<int32_t>> get_subgraphs_rwd(const std::vector<RefPoint>& local_bounds, const int k, const int time_limit) {
    if (local_bounds.size() == 0) {
      return {};
    }
    if (k == 1) {
      std::vector<std::vector<int32_t>> subgraphs;
      for (size_t i = 0; i < local_bounds.size(); ++i) {
        subgraphs.push_back({static_cast<int32_t>(i)});
      }
      return subgraphs;
    }
    auto edges = get_graph_edges(local_bounds);
    auto V = (int)local_bounds.size();
    auto E = (int)edges.size();
    return rwd(V, E, edges, k, time_limit);
  }
};

#endif  // SUBGRAPHS_HPP
