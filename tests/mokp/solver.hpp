#ifndef GHPS_EXPERIMENT_SOLVER_HPP
#define GHPS_EXPERIMENT_SOLVER_HPP

#include <exception>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <ghps/core/timer.hpp>
#include <ghps/mokp/instance.hpp>
#include <ghps/mokp/solvers/branch_bound.hpp>
#include <ghps/mokp/solvers/log.hpp>
#include <ghps/mokp/solvers/milp.hpp>
#include <ghps/mokp/solvers/miqp.hpp>

#include "config.hpp"
#include "helper.hpp"
#include "refpoints.hpp"

namespace moco::ghps {

struct ExperimentRunner {
  const MOKP& instance;
  const TestConfig& config;
  const hv_t nset_hv;
  std::map<Solver, std::ofstream> csv_files;
  std::exception_ptr failure;

  ExperimentRunner(const MOKP& instance, const TestConfig& config)
      : instance(instance), config(config), nset_hv(instance.get_nondominated_set_hv()) {
    for (const auto& [solver, enabled] : config.solvers) {
      if (enabled) init_csv_file(csv_files[solver], solver_csv_filenames.at(solver), csv_header);
    }
  }

  MOKPSolution solve_instance(Solver solver, const std::vector<RefPoint>& refs,
                             const std::vector<std::vector<int64_t>>& coordinates,
                             const std::vector<int64_t>& union_point) const {
    switch (solver) {
      case Solver::Milp:
        return MOKPMILP(instance, {.time_limit = 600.0, .threads = config.solver_threads}).solve(refs);
      case Solver::Miqp:
        return MOKPMIQP(instance, coordinates,
            {.time_limit = 600.0, .threads = config.solver_threads}).solve();
      case Solver::Log:
        return MOKPLogExp(instance.items, instance.N, instance.M, instance.W, coordinates,
            {.time_limit = 600.0, .threads = config.solver_threads}).solve();
      case Solver::BranchBound:
        return MOKPBranchBound::branch_and_bound_greedy(instance,
            RefPoint("up", union_point), refs, Timer(600.0)).solution;
      case Solver::Greedy: {
        auto solution = MOKPBranchBound::compute_greedy_upper_bound(instance, refs);
        MOKPSolution result;
        result.values = std::move(solution.coordinates);
        result.hypervolume = solution.hypervolume;
        result.status = SolverStatus::Optimal;
        return result;
      }
    }
    throw std::invalid_argument("Unknown experiment solver");
  }

  void run(const std::vector<RefPoint>& refs) {
    const Solution simulated = MOKPBranchBound::simulate_solution(instance, refs, {});
    const auto union_point = get_union_point(refs);
    const auto coordinates = ref_to_v(refs);
    const auto refs_str = v_to_str(refs);
    const auto diff_str = diff_v(union_point, simulated.coordinates);
    const auto simulated_str = v_to_str(simulated.coordinates);
    for (auto& [solver, csv_file] : csv_files) {
      try {
        Timer time;
        const auto result = solve_instance(solver, refs, coordinates, union_point);
        const double elapsed = time.elapsed();
        SolverStatus status = result.status;
        std::string message = result.message;
        const bool completed = status == SolverStatus::Optimal || status == SolverStatus::Infeasible;
        if (solver != Solver::Greedy && completed && result.hypervolume != simulated.hypervolume) {
          status = SolverStatus::Error;
          message = "Scalarization hypervolume disagreement: reported=" + moco::ghps::to_string(result.hypervolume)
              + ", simulated=" + moco::ghps::to_string(simulated.hypervolume);
        }
        write_line_to_csv(csv_file, instance.M, instance.N, config.seed, refs.size(),
            refs_str, nset_hv, diff_str, simulated_str, simulated.hypervolume,
            v_to_str(result.values), result.hypervolume, result.hypervolume,
            elapsed, result.nodes, status, result.backend_status, message);
        if (!completed || status == SolverStatus::Error)
          throw std::runtime_error(std::string("Scalarization ended with ") + status_name(status) + ": " + message);
      } catch (...) {
        if (!failure) failure = std::current_exception();
      }
    }
  }

  void finish() const {
    if (failure) std::rethrow_exception(failure);
  }
};

} // namespace moco::ghps

#endif // GHPS_EXPERIMENT_SOLVER_HPP
