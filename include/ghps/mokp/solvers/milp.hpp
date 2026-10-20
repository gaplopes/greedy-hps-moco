#ifndef MOKPMILP_HPP
#define MOKPMILP_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <functional>
#include <set>
#include <tuple>
#include <ghps/mokp/instance.hpp>
#include <ghps/mokp/solvers/branch_bound.hpp>
#include <ghps/mokp/solvers/miqp.hpp>

#include "gurobi_c++.h"
#include <ghps/core/timer.hpp>

namespace moco::ghps {


struct MOKPMILPOptions {
  double time_limit = 600.0;
  int threads = 1;
  Timer deadline = {};
};

class MOKPMILP {
 public:
  MOKPMILP(const MOKP& problem, MOKPMILPOptions options)
      : budget(options.time_limit, options.deadline), items(problem.items),
        N(problem.N),
        M(problem.M),
        W(problem.W),
        time_limit(options.time_limit),
        threads(options.threads), upper_bound(problem.dantzig_upper_bound()) {
    if (!problem.is_maximization) throw std::invalid_argument("MILP supports maximization only");
  }

  static hv_t estimated_variables(int N, int M) {
    hv_t count = 0, combination = 1;
    for (int d = 1; d <= std::min(N, M); ++d) {
      hv_t product;
      if (__builtin_mul_overflow(combination, hv_t(N - d + 1), &product))
        return std::numeric_limits<hv_t>::max();
      combination = product / d;
      if (__builtin_add_overflow(count, combination, &count)) return std::numeric_limits<hv_t>::max();
    }
    return count;
  }

  MOKPSolution solve(const std::vector<RefPoint>& refs) {
    MOKPSolution result;
    try {
      budget.check();
      std::vector<std::vector<int64_t>> coordinates;
      for (const auto& ref : refs) coordinates.push_back(ref.coordinates);
      const hv_t size = estimated_variables(N, M);
      // shortcut: cap polynomial MILP at 2,000 variables; raise only after measuring construction and solve cost.
      if (size > 2000) {
        result = MOKPMIQP(items, N, M, W, coordinates,
            {.time_limit = time_limit, .norm_vec = MOKPMIQP::precision_optimal_norm(upper_bound, 1e-6, budget),
             .threads = threads, .deadline = budget}).solve();
        result.message = "MILP requires " + moco::ghps::to_string(size) + " variables; used MIQP" +
            (result.message.empty() ? "" : ": " + result.message);
        return result;
      }
      GRBEnv env(true);
      env.set(GRB_IntParam_OutputFlag, 0);
      env.set(GRB_IntParam_LogToConsole, 0);
      env.set(GRB_IntParam_Threads, threads);
      env.start();
      GRBModel model(env);
      model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
      model.set(GRB_DoubleParam_TimeLimit, time_limit);
      model.set(GRB_DoubleParam_MIPGap, 1e-6);
      model.set(GRB_IntParam_NumericFocus, 3);
      solveGeneral(model, refs, result);
    } catch (GRBException& e) {
      result.status = SolverStatus::Error;
      result.message = e.getMessage();
    } catch (const TimeLimitReached& e) {
      result.status = SolverStatus::TimeLimit;
      result.message = e.what();
    } catch (const std::exception& e) {
      result.status = SolverStatus::Error;
      result.message = e.what();
    }
    return result;
  }

 private:
  Timer budget;
  // Problem data
  std::vector<Item> items;
  int32_t N;  // items
  int32_t M;  // objectives
  int64_t W;  // capacity
  // Gurobi parameters
  double time_limit;
  int threads;
  std::vector<double> upper_bound;
  std::vector<std::vector<int64_t>> refs;

  void solveGeneral(GRBModel& model, const std::vector<RefPoint>& refPoints, MOKPSolution& result) {

    const int M = this->M;
    const int N = this->N;
    const int64_t MIN = std::numeric_limits<int64_t>::min();

    std::vector<int64_t> allMax(M, MIN);
    for (const auto& rp : refPoints) {
      budget.check();
      for (int d = 0; d < M; ++d) allMax[d] = std::max(allMax[d], rp.coordinates[d]);
    }

    const std::vector<std::vector<int64_t>> refs = [&]() {
      std::vector<std::vector<int64_t>> r(refPoints.size(), std::vector<int64_t>(M));
      for (size_t i = 0; i < refPoints.size(); ++i) {
        budget.check();
        for (int d = 0; d < M; ++d) r[i][d] = refPoints[i].coordinates[d];
      }
      return r;
    }();

    std::vector<std::vector<int64_t>> C(M, std::vector<int64_t>(N));
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        C[i][j] = items[j].values[i];
    
    // ── 1. Precompute inclusion-exclusion subsets ─────────────────────────
    const int K = static_cast<int>(refs.size());
    struct SubsetMax {
      std::vector<int64_t> r;
      int subsetSize;
    };
    std::vector<SubsetMax> subsets;
    for (int mask = 1; mask < (1 << K); ++mask) {
      budget.check();
      std::vector<int64_t> r(M, std::numeric_limits<int64_t>::min());
      int bits = 0;
      for (int m = 0; m < K; ++m) {
        budget.check();
        if (mask & (1 << m)) {
          ++bits;
          for (int d = 0; d < M; ++d) r[d] = std::max(r[d], refs[m][d]);
        }
      }
      subsets.push_back({std::move(r), bits});
    }

    // ── 2. Combination generators ────────────────────────────────────────
    auto for_each_combination = [&](int n, int m, bool repeated, const auto& visit) {
      std::vector<int> current;
      std::function<void(int)> next = [&](int start) {
        budget.check();
        if (static_cast<int>(current.size()) == m) { visit(current); return; }
        const int end = repeated ? n : n - (m - static_cast<int>(current.size())) + 1;
        for (int i = start; i < end; ++i) {
          budget.check();
          current.push_back(i);
          next(repeated ? i : i + 1);
          current.pop_back();
        }
      };
      next(0);
    };

    // ── 3. Universal polynomial coefficient generator ────────────────────
    //
    // calc_coeff(t, r) computes the coefficient of the monomial x_{t_0}·…·x_{t_{m-1}}
    // arising from expanding Π_d (f_d(x) − r_d).
    //
    // Algorithm:
    //   For each way to choose m dimensions out of M (bitmask over dims):
    //     - The chosen dimensions take variable factors from the tuple.
    //     - The remaining (M − m) dimensions contribute r_d factors.
    //     - Sum over all permutations of the tuple across chosen dimensions.
    //
    // Example (M=3, tuple = {i,j}, r = {r0,r1,r2}):
    //   Bitmask 011 → dims {0,1} for vars, dim {2} for r:
    //     r2 · Σ_{perm(i,j)} C[0][π0]·C[1][π1]
    //     = r2 · (C[0][i]·C[1][j] + C[0][j]·C[1][i])
    //   Bitmask 101 → dims {0,2}, dim {1}: r1 · (...)
    //   Bitmask 110 → dims {1,2}, dim {0}: r0 · (...)
    //
    auto calc_coeff = [&](const std::vector<int>& var_idx,
                          const std::vector<int64_t>& r) -> hv_t {
      hv_t total = 0;
      int num_vars = static_cast<int>(var_idx.size());

      for (int mask = 0; mask < (1 << M); ++mask) {
        budget.check();
        int bits = 0;
        for (int d = 0; d < M; ++d)
          if (mask & (1 << d)) bits++;
        if (bits != num_vars) continue;

        std::vector<int> var_dims, r_dims;
        for (int d = 0; d < M; ++d) {
          budget.check();
          if (mask & (1 << d)) var_dims.push_back(d);
          else                 r_dims.push_back(d);
        }

        hv_t r_prod = 1;
        for (int d : r_dims) r_prod *= r[d];

        std::vector<int> perm = var_idx;
        std::sort(perm.begin(), perm.end());
        do {
          budget.check();
          hv_t ip = 1;
          for (size_t p = 0; p < perm.size(); ++p)
            ip *= C[var_dims[p]][perm[p]];
          total += ip * r_prod;
        } while (std::next_permutation(perm.begin(), perm.end()));
      }
      return total;
    };

    // ── 4. Variable generation (strict index tuples only) ────────────────
    //
    // Since x_j ∈ {0,1} ⇒ x_j² = x_j, any monomial with repeated indices
    // collapses to a lower-degree monomial over the *unique* indices.
    // Therefore we only create linearization variables for strictly
    // increasing tuples {i < j}, {i < j < k}, etc.
    //
    // The map key is the sorted tuple; degree-1 maps to the original x_j.
    //
    std::map<std::vector<int>, GRBVar> vars;

    // Degree-1 variables are the decision variables themselves
    for (int i = 0; i < N; ++i) {
      budget.check();
      vars[{i}] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + std::to_string(i));
    }

    // Higher-degree linearization variables: y_{ij}, w_{ijk}, ...
    for (int m = 2; m <= M; ++m) {
      budget.check();
      for_each_combination(N, m, false, [&](const auto& t) {
        vars[t] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "lin_" + tuple_name(t));
      });
    }

    // ── 5. Objective construction ────────────────────────────────────────
    //
    // Assemble:  obj = Σ_{m=0}^{M} (−1)^{M−m} · Σ_t coeff_sum(t) · vars[collapse(t)]
    //
    // For the highest-degree terms (m = M), the reference part is empty,
    // so coeff is independent of the reference points. For lower degrees
    // (m < M), inclusion–exclusion over subsets modulates the coefficients.
    //
    // The "collapse" step maps tuples with repeats to their unique-index
    // equivalent (e.g., {0,0,1} → {0,1}), exploiting x² = x for binaries.
    //
    GRBLinExpr obj = 0;
    std::vector<int64_t> zero_r(M, 0);

    for (int m = 1; m <= M; ++m) {
      budget.check();
      int sign_M = ((M - m) % 2 == 0) ? +1 : -1;  // (−1)^{M−m}
      for_each_combination(N, m, true, [&](const auto& t) {
        hv_t coeff_sum = 0;
        if (m == M) {
          coeff_sum = calc_coeff(t, zero_r);
        } else {
          for (const auto& sub : subsets) {
            budget.check();
            int sign_sub = (sub.subsetSize % 2 == 1) ? +1 : -1;
            coeff_sum += sign_sub * calc_coeff(t, sub.r);
          }
        }
        coeff_sum *= sign_M;

        if (coeff_sum != 0) {
          // Collapse repeated vars: x_0·x_0·x_1 → y_{0,1} (since x_i²=x_i for binary)
          std::set<int> unique(t.begin(), t.end());
          std::vector<int> collapsed(unique.begin(), unique.end());
          obj += static_cast<double>(coeff_sum) * vars[collapsed];
        }
      });
    }

    // Constant term (m = 0)
    hv_t const_sum = 0;
    int sign_M_0 = (M % 2 == 0) ? +1 : -1;
    for (const auto& sub : subsets) {
      budget.check();
      int sign_sub = (sub.subsetSize % 2 == 1) ? +1 : -1;
      const_sum += sign_sub * calc_coeff({}, sub.r);
    }
    obj += static_cast<double>(sign_M_0 * const_sum);

    model.setObjective(obj, GRB_MAXIMIZE);

    // ── 6. McCormick envelopes ───────────────────────────────────────────
    //
    // For each linearization variable (representing a product of binaries),
    // add McCormick envelope constraints that make the LP relaxation tight.
    //
    // General rule for z ≈ a · b where a, b ∈ [0,1]:
    //   z ≤ a,  z ≤ b,  z ≥ a + b − 1,  z ≥ 0
    //
    // These are exact (not a relaxation) when a, b are binary.
    //
    // For degree-m variables, we decompose into (degree-(m-1)) × (degree-1)
    // and apply McCormick to each decomposition, giving 3·m constraints
    // per variable. Multiple decompositions tighten the LP relaxation.
    //
    for (int m = 2; m <= M; ++m) {
      budget.check();
      for_each_combination(N, m, false, [&](const auto& t) {
        if (m == 2) {
          // For degree-2, the two possible splits {i}×{j} and {j}×{i}
          // produce identical McCormick constraints (product is symmetric).
          // Emit each constraint exactly once:
          //   lin_{ij} <= x_i,  lin_{ij} <= x_j,
          //   lin_{ij} >= x_i + x_j - 1
          std::string tag = "mc_" + tuple_name(t);
          model.addConstr(vars[t] <= vars[{t[0]}], tag + "_ub1");
          model.addConstr(vars[t] <= vars[{t[1]}], tag + "_ub2");
          model.addConstr(vars[t] >= vars[{t[0]}] + vars[{t[1]}] - 1,
                        tag + "_lb");
        } else {
          // For degree >= 3, each split through a different sub-tuple
          // produces distinct constraints, so emit all of them.
          for (size_t p = 0; p < t.size(); ++p) {
            budget.check();
            std::vector<int> t_minus_p = t;
            t_minus_p.erase(t_minus_p.begin() + static_cast<long>(p));
            int x_idx = t[p];

            std::string tag = "mc_" + tuple_name(t) + "_" + std::to_string(p);
            model.addConstr(vars[t] <= vars[t_minus_p], tag + "_ub1");
            model.addConstr(vars[t] <= vars[{x_idx}],   tag + "_ub2");
            model.addConstr(vars[t] >= vars[t_minus_p] + vars[{x_idx}] - 1,
                          tag + "_lb");
          }
        }
      });
    }

    // ── 6.5. Capacity and Reference Point Constraints ────────────────────
    GRBLinExpr cap = 0;
    for (int i = 0; i < N; ++i) cap += items[i].weight * vars[{i}];
    model.addConstr(cap <= W, "Constraint_W");

    std::vector<GRBLinExpr> f(M, 0);
    for (int i = 0; i < N; ++i) {
      budget.check();
      for (int d = 0; d < M; ++d) f[d] += items[i].values[d] * vars[{i}];
    }
    for (int d = 0; d < M; ++d) {
      budget.check();
      model.addConstr(f[d] >= static_cast<double>(hv_t(allMax[d]) + 1));
    }

    // --- 7. Optimize & Extract ---
    model.update();
    budget.check();
    model.set(GRB_DoubleParam_TimeLimit, budget.left());
    model.optimize();

    result.backend_status = model.get(GRB_IntAttr_Status);
    switch (result.backend_status) {
      case GRB_OPTIMAL: result.status = SolverStatus::Optimal; break;
      case GRB_INFEASIBLE: result.status = SolverStatus::Infeasible; break;
      case GRB_TIME_LIMIT: result.status = SolverStatus::TimeLimit; break;
      case GRB_INTERRUPTED: result.status = SolverStatus::Interrupted; break;
      case GRB_NUMERIC: result.status = SolverStatus::NumericalFailure; break;
      case GRB_UNBOUNDED: result.status = SolverStatus::Unbounded; break;
      default: result.status = SolverStatus::Other; break;
    }
    result.nodes = static_cast<int64_t>(model.get(GRB_DoubleAttr_NodeCount));
    if (model.get(GRB_IntAttr_SolCount) == 0) return;
    result.items.assign(N, 0);
    result.values.assign(M, 0);
    for (int i = 0; i < N; ++i) {
      result.items[items[i].idx] = vars[{i}].get(GRB_DoubleAttr_X) > 0.5;
      if (!result.items[items[i].idx]) continue;
      result.weight += items[i].weight;
      for (int d = 0; d < M; ++d) result.values[d] += items[i].values[d];
    }
    result.hypervolume = static_cast<hv_t>(std::round(model.get(GRB_DoubleAttr_ObjVal)));
  }

  static std::string tuple_name(const std::vector<int>& t) {
    std::string s;
    for (size_t i = 0; i < t.size(); ++i) {
      if (i > 0) s += "_";
      s += std::to_string(t[i]);
    }
    return s;
  }

};

} // namespace moco::ghps

#endif  // MOKPMILP_HPP
