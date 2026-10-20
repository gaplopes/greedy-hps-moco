#ifndef GHPS_EXPERIMENT_CONFIG_HPP
#define GHPS_EXPERIMENT_CONFIG_HPP

#include <cstdint>
#include <map>

namespace moco::ghps {

enum class Solver { Milp, Miqp, Log, Greedy, BranchBound };

struct TestConfig {
  int32_t seed = 42;
  std::map<Solver, bool> solvers = {
      {Solver::Milp, false}, {Solver::Miqp, true}, {Solver::Log, false},
      {Solver::Greedy, false}, {Solver::BranchBound, false}};
  int solver_threads = 4;
  int num_ref_points = 30;
  int max_k = 5;
  int max_permutations = 30;
};

} // namespace moco::ghps

#endif // GHPS_EXPERIMENT_CONFIG_HPP
