#pragma once

// Reverse search with a dictionary: connected subsets differ by one vertex swap.
// Komusiewicz and Sommer (SOFSEM 2019), Section 1:
// https://www.uni-marburg.de/de/fb12/arbeitsgruppen/algorith/paper/enucon.pdf/@@download/file/enucon.pdf

#include <iostream>
#include <map>
#include <utility>
#include <cassert>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <algorithm>
#include <ghps/core/timer.hpp>

namespace moco::ghps {

struct NodeAux {
  std::string id;
  std::set<int> neighbors_index;

  NodeAux(std::string id) : id(id) {}

  bool operator<(const NodeAux &other) const {
    return id < other.id;
  }

  bool operator>(const NodeAux &other) const {
    return id > other.id;
  }

  bool operator==(const NodeAux &other) const {
    return id == other.id;
  }

  void add_neighbor(int index) {
    neighbors_index.insert(index);
  }

  bool is_neighbor(int index) {
    return neighbors_index.find(index) != neighbors_index.end();
  }

  std::string to_string() const {
    std::string str = id + " {";
    for (auto it = neighbors_index.begin(); it != neighbors_index.end(); ++it) {
      if (it != neighbors_index.begin()) {
        str += ", ";
      }
      str += std::to_string(*it);
    }
    str += "}";
    return str;
  }
};

class UnionFind {
 public:
  UnionFind(const std::vector<int> &nodes) {
    int n = nodes.size();
    parent.resize(n);
    rank.resize(n, 0);
    for (int i = 0; i < n; ++i) {
      parent[i] = i;
      name_to_index[nodes[i]] = i;
    }
    name = nodes;
  }

  int findIndex(const int &name) {
    return name_to_index[name];
  }

  int find(int index) {
    if (index != parent[index]) {
      parent[index] = find(parent[index]);
    }
    return parent[index];
  }

  void unionSets(const int ind_a, const int ind_b) {
    int idx_a = find(findIndex(ind_a));
    int idx_b = find(findIndex(ind_b));

    if (idx_a == idx_b) {
      return;
    }

    if (rank[idx_a] > rank[idx_b]) {
      parent[idx_b] = idx_a;
    } else {
      parent[idx_a] = idx_b;
      if (rank[idx_a] == rank[idx_b]) {
        rank[idx_b]++;
      }
    }
  }

  void printParent() {
    std::cout << "parent: ";
    for (int i : parent) {
      std::cout << i << " ";
    }
    std::cout << std::endl;

    std::cout << "names: ";
    for (auto &entry : name_to_index) {
      std::cout << entry.first << ":" << entry.second << " ";
    }
    std::cout << std::endl;
  }

  std::vector<int> representatives() {
    std::vector<int> representatives;
    for (int i = 0; i < (int)name.size(); ++i) {
      if (i == parent[i]) {
        representatives.push_back(name[i]);
      }
    }
    return representatives;
  }

  std::vector<std::unordered_set<int>> components() {
    std::vector<int> reps = representatives();

    for (int i = 0; i < (int)reps.size(); ++i) {
      reps[i] = findIndex(reps[i]);
    }

    std::vector<std::unordered_set<int>> compos(reps.size());

    for (int i = 0; i < (int)name.size(); ++i) {
      auto it = std::lower_bound(reps.begin(), reps.end(), find(i));
      compos[it - reps.begin()].insert(name[i]);
    }

    return compos;
  }

 private:
  std::vector<int> parent;
  std::vector<int> rank;
  std::unordered_map<int, int> name_to_index;
  std::vector<int> name;
};

struct Graph {
  std::vector<NodeAux> nodes;
  std::set<std::string> names;
  std::vector<std::string> vertices_map;
  std::vector<std::pair<int, int>> edges;
  std::unordered_map<std::string, int> id_to_index;

  Graph(std::vector<NodeAux> nodes,
        std::set<std::string> names,
        std::vector<std::string> vertices_map,
        std::vector<std::pair<int, int>> edges,
        std::unordered_map<std::string, int> id_to_index) : nodes(std::move(nodes)),
                                                            names(std::move(names)),
                                                            vertices_map(std::move(vertices_map)),
                                                            edges(std::move(edges)),
                                                            id_to_index(std::move(id_to_index)) {}

  static auto from_input(int V, int E, const std::vector<std::pair<int, int>> &edges, const Timer& deadline = Timer()) -> Graph {
    std::vector<NodeAux> nodes;
    std::vector<std::string> vertices_map;
    for (int i = 0; i < V; i++) {
      deadline.check();
      nodes.emplace_back(std::to_string(i));
      vertices_map.push_back(std::to_string(i));
    }

    std::set<std::string> names;
    std::vector<std::pair<int, int>> edges_aux;
    std::unordered_map<std::string, int> id_to_index;

    for (int i = 0; i < E; i++) {
      deadline.check();
      std::string u = std::to_string(edges[i].first);
      std::string v = std::to_string(edges[i].second);
      id_to_index[u] = edges[i].first;
      id_to_index[v] = edges[i].second;
      nodes[id_to_index[u]].add_neighbor(id_to_index[v]);
      nodes[id_to_index[v]].add_neighbor(id_to_index[u]);
      edges_aux.emplace_back(id_to_index[u], id_to_index[v]);
    }

    return Graph(nodes, names, vertices_map, edges_aux, id_to_index);
  }

  static auto from_search_tree(std::vector<NodeAux> nodes,
                               std::set<std::string> names,
                               std::vector<std::string> vertices_map,
                               const std::vector<std::pair<int, int>> &edges,
                               std::unordered_map<std::string, int> id_to_index, const Timer& deadline = Timer()) -> Graph {
    auto old_nodes = nodes;
    nodes.clear();
    for (int i = 0; i < (int)vertices_map.size(); i++) {
      deadline.check();
      nodes.push_back(NodeAux(std::to_string(i)));
    }
    for (int i = 0; i < (int)edges.size(); i++) {
      deadline.check();
      auto edge = edges[i];
      int u = edge.first;
      int v = edge.second;
      nodes[std::stoi(old_nodes[u].id)].add_neighbor(std::stoi(old_nodes[v].id));
      nodes[std::stoi(old_nodes[v].id)].add_neighbor(std::stoi(old_nodes[u].id));
    }
    return Graph(nodes, names, vertices_map, edges, id_to_index);
  }

  void print_nodes() {
    std::cout << "Nodes: ";
    std::cout << nodes.size() << std::endl;
    for (auto const &n : nodes) {
      std::cout << n.to_string() << std::endl;
    }
  }

  void print_edges() {
    std::cout << "Edges: ";
    std::cout << edges.size() << std::endl;
    for (auto const &e : edges) {
      std::cout << e.first << " " << e.second << std::endl;
    }
  }

  int get_index(std::string id) {
    return id_to_index[id];
  }

  bool are_connected(int index1, int index2) {
    return nodes[index1].is_neighbor(index2);
  }

  std::string ids_to_subgraph(const std::vector<int>& index) {
    std::vector<int> indexes;
    for (auto const &id : index) {
      indexes.push_back(std::stoi(vertices_map[id]));
    }
    std::sort(indexes.begin(), indexes.end());
    std::string str = std::to_string(indexes[0]);
    for (int i = 1; i < (int)indexes.size(); i++) {
      str += " " + std::to_string(indexes[i]);
    }
    return str;
  }

  std::vector<std::vector<NodeAux>> ids_to_nodes(const std::vector<std::string>& subgraphs_k, const Timer& deadline = Timer()) {
    std::vector<std::vector<NodeAux>> subgraphs;
    for (auto const &subgraph_str : subgraphs_k) {
      deadline.check();
      std::vector<NodeAux> subgraph;
      std::stringstream ss(subgraph_str);
      std::string id;
      while (ss >> id) {
        deadline.check();
        int index = std::stoi(id);
        subgraph.push_back(nodes[index]);
      }
      subgraphs.push_back(subgraph);
    }
    return subgraphs;
  }

  std::vector<std::vector<NodeAux>> rwd(const std::vector<std::vector<int>>& components, int k, const Timer& deadline) {
    std::vector<int> vertices_map_inverse(vertices_map.size());
    for (int i = 0; i < (int)vertices_map.size(); i++) {
      deadline.check();
      vertices_map_inverse[std::stoi(vertices_map[i])] = i;
    }
    std::vector<std::string> R;                   // Set of subgraphs of size k
    int reduced_size = k - 1;                     // Size of reduced subgraph
    std::set<std::string> subgraphs_k;            // Set of subgraphs of size k, represented as strings of node ids
    std::queue<std::vector<int>> subgraph_queue;  // Queue of subgraphs of size k
    for (auto const &component : components) {
      deadline.check();
      std::string subgraph_str = ids_to_subgraph(component);  // Convert subgraph to string
      subgraphs_k.insert(subgraph_str);                       // Add subgraph to set
      subgraph_queue.push(component);                         // Add subgraph to queue
      while (!subgraph_queue.empty()) {
        deadline.check();
        std::vector<int> subgraph = subgraph_queue.front();  // Get subgraph from queue
        subgraph_queue.pop();
        R.push_back(ids_to_subgraph(subgraph));  // Add subgraph to set of subgraphs of size k
        std::vector<int> subgraph_copy = subgraph;
        for (auto const &vertex : subgraph_copy) {
          deadline.check();
          subgraph.erase(std::remove_if(subgraph.begin(), subgraph.end(), [&vertex](int const &n) { return n == vertex; }));
          UnionFind UF(subgraph);
          for (int i = 0; i < reduced_size; i++) {
            deadline.check();
            int vertex1 = subgraph[i];
            int vertex1_id = std::stoi(vertices_map[subgraph[i]]);
            for (int j = i + 1; j < reduced_size; j++) {
              deadline.check();
              int vertex2 = subgraph[j];
              int vertex2_id = std::stoi(vertices_map[subgraph[j]]);
              if (are_connected(vertex1_id, vertex2_id)) {
                UF.unionSets(vertex1, vertex2);
              }
            }
          }
          std::vector<std::unordered_set<int>> connected_components = UF.components();
          std::set<int> possible_neighbors;
          for (int i = 0; i < (int)subgraph.size(); i++) {
            deadline.check();
            auto v = subgraph[i];
            auto v_id = std::stoi(vertices_map[v]);
            auto ids_set = std::vector<std::pair<int, int>>();
            for (auto const &neighbor_id : nodes[v_id].neighbors_index) {
              deadline.check();
              ids_set.emplace_back(vertices_map_inverse[neighbor_id], neighbor_id);
            }
            std::sort(ids_set.begin(), ids_set.end());
            for (auto const &neig : ids_set) {
              deadline.check();
              auto neighbor_id = neig.second;
              auto neighbor = neig.first;
              auto in_set = possible_neighbors.find(neighbor);
              auto in_S_copy = std::find_if(subgraph_copy.begin(), subgraph_copy.end(), [&neighbor](int const &n) { return n == neighbor; });
              if (in_set == possible_neighbors.end() && in_S_copy == subgraph_copy.end()) {
                possible_neighbors.insert(neighbor);
                bool is_connected = true;
                for (auto const &component : connected_components) {
                  deadline.check();
                  bool check = false;
                  for (auto const &n : component) {
                    deadline.check();
                    auto n_id = std::stoi(vertices_map[n]);
                    bool is_neighbor = are_connected(n_id, neighbor_id);
                    if (is_neighbor) {
                      check = true;
                      break;
                    }
                  }
                  if (!check) {
                    is_connected = false;
                    break;
                  }
                }
                if (!is_connected) {
                  continue;
                }
                subgraph.push_back(neighbor);
                std::string S_str = ids_to_subgraph(subgraph);
                if (subgraphs_k.find(S_str) == subgraphs_k.end()) {
                  subgraphs_k.insert(S_str);
                  subgraph_queue.push(subgraph);
                }
                subgraph.pop_back();
              }
            }
          }
          subgraph.push_back(vertex);
        }
      }
    }
    return ids_to_nodes(R, deadline);
  }
};

struct SearchTree {
  Graph graph;                               // Vector of nodes
  std::vector<std::vector<int>> components;  // Vector of components

  SearchTree(Graph graph) : graph(std::move(graph)), components(std::vector<std::vector<int>>()) {}

  Graph create_search_tree(int k, const Timer& deadline = Timer()) {



    std::map<std::int64_t, std::int64_t> dict_ids;         // Map of subgraph ids to number of occurrences
    std::vector<std::string> vertex_names;                 // Vector of vertex names
    std::vector<std::vector<int32_t>> components_indexes;  // Vector of components
    int32_t n = graph.nodes.size();                        // Number of nodes
    int32_t id = n - 1;                                    // Current subgraph id
    std::vector<bool> visited(n, false);                   // Vector of visited nodes
    for (int i = 0; i < n; i++) {
      deadline.check();
      if (visited[i]) continue;
      std::set<int32_t> component;                             // Set of nodes in current component
      std::stack<int32_t> s;                                   // Queue of nodes to visit
      s.push(i);                                               // Add current node to queue
      int32_t smallest_id = n - (int32_t)vertex_names.size();  // Smallest id in current component
      while (!s.empty()) {
        deadline.check();
        int32_t u = s.top();  // Get node from queue
        s.pop();              // Remove node from queue
        if (component.find(u) != component.end()) continue;     // If node is already in component, skip
        component.insert(u);                                    // Add node to component
        visited[u] = true;                                      // Mark node as visited
        dict_ids[u] = id;                                       // Add node to subgraph id
        for (auto const &v : graph.nodes[u].neighbors_index) {  // Add neighbors to queue
          deadline.check();
          if (component.find(v) == component.end()) {
            s.push(v);
          }
        }

        vertex_names.insert(vertex_names.begin(), graph.nodes[u].id);  // Insert vertex name into first position of vector
        if ((int)component.size() == k) {                              // If component size is k, add component to vector of components
          std::vector<int> component_indexes;                          // The smallest id - 1 is the id of the last node added to the component
          for (int j = smallest_id - k; j < smallest_id; j++) {
            deadline.check();
            component_indexes.push_back(j);
          }
          components_indexes.push_back(component_indexes);
        }
      }
    }

    std::vector<std::pair<int32_t, int32_t>> edges;
    for (auto &edge : graph.edges) {
      deadline.check();
      std::string source = graph.vertices_map[edge.first];
      std::string target = graph.vertices_map[edge.second];
      int32_t source_index = find(vertex_names.begin(), vertex_names.end(), source) - vertex_names.begin();
      int32_t target_index = find(vertex_names.begin(), vertex_names.end(), target) - vertex_names.begin();
      int32_t new_source = graph.id_to_index[vertex_names[target_index]];
      int32_t new_target = graph.id_to_index[vertex_names[source_index]];
      edges.push_back(std::make_pair(new_source, new_target));
    }

    Graph final_graph = Graph::from_search_tree(graph.nodes, graph.names, vertex_names, edges, graph.id_to_index, deadline);
    components = components_indexes;
    return final_graph;
  }
};

inline std::vector<std::vector<int32_t>> rwd(const int V, const int E, const std::vector<std::pair<int32_t, int32_t>> &edges, int k, const Timer& deadline) {
  deadline.check();
  auto G = Graph::from_input(V, E, edges, deadline);
  // In this polynomial delay algorithm the vertices have to be sorted according to a DFS search tree
  auto search_tree = SearchTree(std::move(G));
  G = search_tree.create_search_tree(k, deadline);
  if (search_tree.components.empty()) return {};
  auto R = G.rwd(search_tree.components, k, deadline);
  std::vector<std::vector<int32_t>> R_int;
  for (auto const &subgraph : R) {
    deadline.check();
    std::vector<int32_t> subgraph_int;
    for (auto const &node : subgraph) {
      deadline.check();
      subgraph_int.push_back(stoi(node.id));
    }
    R_int.push_back(subgraph_int);
  }
  return R_int;
}

inline std::vector<std::vector<int32_t>> rwd(int V, int E, const std::vector<std::pair<int32_t, int32_t>>& edges, int k, double time_limit) {
  return rwd(V, E, edges, k, Timer(time_limit));
}

} // namespace moco::ghps
