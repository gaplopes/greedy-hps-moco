#ifndef LOG_HPP
#define LOG_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ghps/mokp/instance.hpp>
#include <memory>
#include <numeric>
#include <vector>

#include "gurobi_c++.h"
#include <ghps/core/timer.hpp>

namespace moco::ghps {


struct MOKPLogExpOptions {
  double time_limit = 600.0;
  int threads = 1;
  Timer deadline = {};
};

/**
 * Unified MOKP model structure that can handle both single and multiple reference points
 */
struct MOKPLogExp {
 private:
  Timer budget;
  MOKPSolution initial_result;
  std::unique_ptr<GRBEnv> env;
  std::unique_ptr<GRBModel> model;
  std::vector<Item> items;         // Items from the MOKP instance
  int32_t N;                       // Number of items
  int32_t M;                       // Number of objectives
  int64_t W;                       // Knapsack capacity
  std::vector<GRBVar> x;           // Decision variables for items
  std::vector<GRBLinExpr> f_expr;  // Linear expressions for objectives

  // Common variables
  GRBConstr capacity_constr;  // Knapsack capacity constraint

  // Reference Point Data
  std::vector<std::vector<int64_t>> refs;  // Reference points for optimization

  // Model parameters
  double time_limit = 600.0;                     // In seconds
  static constexpr double EPS = 1e-12;           // Minimum value for log argument (stability)
  static constexpr double DOMINATION_EPS = 0.5;  // Minimum value for strictly dominating solutions (integers)

  // Returns a GRBVar S (continuous) representing the sum of logs: construct S = sum_i log(f_i(x) - ref[i]) for a given ref
  GRBVar build_sum_of_logs(const std::vector<int64_t>& ref, const std::string& tag) {
    std::vector<GRBVar> z(M), y(M);

    for (int i = 0; i < M; ++i) {
      budget.check();
      // z = f_i(x) - ref[i]
      z[i] = model->addVar(DOMINATION_EPS, GRB_INFINITY, 0.0, GRB_CONTINUOUS, tag + "_z_" + std::to_string(i));
      model->addConstr(z[i] == f_expr[i] - ref[i], tag + "_zdef_" + std::to_string(i));

      // y_i = log(z_i)
      y[i] = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, tag + "_y_" + std::to_string(i));
      model->addGenConstrLog(z[i], y[i], tag + "_log_" + std::to_string(i));
    }

    // S = sum_i y_i
    GRBVar S = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, tag + "_S");
    GRBLinExpr sumy = 0.0;
    for (int i = 0; i < M; ++i) sumy += y[i];
    model->addConstr(S == sumy, tag + "_Sdef");
    return S;
  }

  // Given vars v[0..k-1] compute LSE = log(sum_i exp(v[i])).
  GRBVar log_sum_exp(const std::vector<GRBVar>& vars, const std::string& tag) {
    if (vars.empty()) throw std::runtime_error("log_sum_exp: empty var list");

    // maxVar = max(vars)
    GRBVar maxVar = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, tag + "_max");
    // Gurobi expects pointer to array of GRBVar; pass address of first element
    model->addGenConstrMax(maxVar, const_cast<GRBVar*>(vars.data()), static_cast<int>(vars.size()), -GRB_INFINITY, tag + "_maxConstr");

    // For each var: diff = var - maxVar ; e = exp(diff)
    std::vector<GRBVar> expVars;
    expVars.reserve(vars.size());
    for (size_t k = 0; k < vars.size(); ++k) {
      budget.check();
      GRBVar diff = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, tag + "_diff_" + std::to_string(k));
      model->addConstr(diff == vars[k] - maxVar, tag + "_diff_def_" + std::to_string(k));
      // Since diff <= 0 by construction, exp(diff) ∈ (0, 1]; tighten upper bound to 1.0
      GRBVar e = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, tag + "_exp_" + std::to_string(k));
      model->addGenConstrExp(diff, e, tag + "_expConstr_" + std::to_string(k));
      expVars.push_back(e);
    }

    // sumExp = sum(e)
    // At least one term equals maxVar, so one exp(diff)=1, and others ∈ (0,1];
    // therefore sumExp ∈ [1, vars.size()].
    GRBVar sumExp = model->addVar(1.0, static_cast<double>(vars.size()), 0.0, GRB_CONTINUOUS, tag + "_sumExp");
    GRBLinExpr tmp = 0.0;
    for (auto& ev : expVars) tmp += ev;
    model->addConstr(sumExp == tmp, tag + "_sumExp_def");

    // logSum = log(sumExp)
    GRBVar logSum = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, tag + "_logSum");
    model->addGenConstrLog(sumExp, logSum, tag + "_logSumConstr");

    // result = maxVar + logSum
    GRBVar result = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, tag + "_LSE");
    model->addConstr(result == maxVar + logSum, tag + "_LSE_def");
    return result;
  }

 public:
  MOKPLogExp(
      const std::vector<Item>& items,
      const int32_t& N,
      const int32_t& M,
      const int64_t& W,
      const std::vector<std::vector<int64_t>>& refs,
      MOKPLogExpOptions options)
      : budget(options.time_limit, options.deadline), items(items), N(N), M(M), W(W), refs(refs), time_limit(options.time_limit) {
    // Initialize Gurobi environment and model
    try {
      budget.check();

      env = std::make_unique<GRBEnv>(true);
      env->set(GRB_IntParam_OutputFlag, 0);
      env->set(GRB_IntParam_Threads, options.threads);
      env->start();
      // Create the model from the started environment.
      model = std::make_unique<GRBModel>(*env);
      model->set(GRB_DoubleParam_TimeLimit, time_limit);
      model->set(GRB_DoubleParam_MIPGap, 1e-6);
      model->set(GRB_IntParam_NumericFocus, 3);
      model->set(GRB_IntParam_DualReductions, 0);  // Resolve infeasible status without an input precheck.
      model->set(GRB_IntParam_FuncNonlinear, 1);  // Enable non-linear functions

      // Adaptive FeasibilityTol:
      // k=1: High precision (1e-9) is needed for the steep log gradients.
      // k>1: Standard precision (1e-6) prevents noise accumulation in inclusion-exclusion.
      if (refs.size() == 1) {
        model->set(GRB_DoubleParam_FeasibilityTol, 1e-9);
      } else {
        model->set(GRB_DoubleParam_FeasibilityTol, 1e-6);
      }



      // Add decision variables x_j for all items
      x.resize(N);
      for (int j = 0; j < N; ++j) {
        budget.check();
        x[j] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + std::to_string(j));
      }

      // Add capacity constraint: sum(w_j * x_j) <= W
      GRBLinExpr capacity_expr = 0;
      for (int j = 0; j < N; ++j) {
        budget.check();
        capacity_expr += items[j].weight * x[j];
      }
      capacity_constr = model->addConstr(capacity_expr <= W, "knapsack_capacity");

      // Precompute linear expressions f_i(x) = sum_j v_ij * x_j
      f_expr.resize(M);
      for (int i = 0; i < M; ++i) {
        budget.check();
        f_expr[i] = 0;
        for (int j = 0; j < N; ++j) {
          budget.check();
          f_expr[i] += items[j].values[i] * x[j];
        }
      }

      int32_t k = static_cast<int32_t>(refs.size());

      switch (k) {
        case 1:
          initialize_model_one_reference_point();
          break;
        case 2:
          initialize_model_two_reference_points();
          break;
        case 3:
          initialize_model_three_reference_points();
          break;
        default:
          initialize_model_general_reference_points();
          break;
      }

    } catch (GRBException& e) {
      std::cerr << "Gurobi Error [" << e.getErrorCode() << "]: "
                << e.getMessage() << std::endl;
      initial_result.status = SolverStatus::Error;
      initial_result.message = e.getMessage();
    } catch (const TimeLimitReached& e) {
      initial_result.status = SolverStatus::TimeLimit;
      initial_result.message = e.what();
    } catch (const std::exception& e) {
      initial_result.status = SolverStatus::Error;
      initial_result.message = e.what();
    }
  }

  void initialize_model_one_reference_point() {
    try {
      // Reference point
      const std::vector<int64_t>& r1 = refs[0];

      // Add auxiliary variables z_i and y_i for each objective
      std::vector<GRBVar> z(M);
      std::vector<GRBVar> y(M);
      std::vector<GRBConstr> z_constrs(M);

      for (int i = 0; i < M; ++i) {
        budget.check();
        z[i] = model->addVar(DOMINATION_EPS, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z_" + std::to_string(i));
        y[i] = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_" + std::to_string(i));
        // Define z_i = f^i(x) - r^i; integer gains must be positive.

        z_constrs[i] = model->addConstr(z[i] == f_expr[i] - r1[i], "def_z_" + std::to_string(i));

        // Add logarithmic constraint: y_i = log(z_i)
        model->addGenConstrLog(z[i], y[i], "log_constr_" + std::to_string(i));
      }

      // Set objective: maximize sum(y_i)
      GRBLinExpr S = 0;
      for (int i = 0; i < M; ++i) S += y[i];
      model->setObjective(S, GRB_MAXIMIZE);
      model->update();
    } catch (GRBException& e) {
      std::cerr << "Gurobi error: " << e.getMessage() << std::endl;
      throw;
    }
  }

  void initialize_model_two_reference_points() {
    try {
      // Reference points
      const std::vector<int64_t>& r1 = refs[0];
      const std::vector<int64_t>& r2 = refs[1];

      std::vector<int64_t> rU(M);
      for (int i = 0; i < M; ++i)
        rU[i] = std::max(r1[i], r2[i]);

      // Build S1, S2, SU
      GRBVar S1 = build_sum_of_logs(r1, "r1");
      GRBVar S2 = build_sum_of_logs(r2, "r2");
      GRBVar SU = build_sum_of_logs(rU, "rU");

      GRBVar A = log_sum_exp({S1, S2}, "pos_lse");
      GRBVar B = SU;

      // Enforce A > B for numerical stability
      model->addConstr(A >= B + EPS, "ensure_A_gt_B");

      // Result = A + log(1 - exp(B - A))
      GRBVar B_minus_A = model->addVar(-GRB_INFINITY, 0.0, 0.0, GRB_CONTINUOUS, "B_minus_A");
      model->addConstr(B_minus_A == B - A, "def_B_minus_A");
      GRBVar exp_B_minus_A = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "exp_B_minus_A");
      model->addGenConstrExp(B_minus_A, exp_B_minus_A, "exp_B_minus_A_constr");
      GRBVar one_minus = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "one_minus");
      model->addConstr(one_minus == 1.0 - exp_B_minus_A, "def_one_minus");
      GRBVar log_one_minus = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "log_one_minus");
      model->addGenConstrLog(one_minus, log_one_minus, "log_one_minus_constr");
      GRBVar result = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "result");
      model->addConstr(result == A + log_one_minus, "def_result");

      // Objective: maximize the stable result
      model->setObjective(GRBLinExpr(result), GRB_MAXIMIZE);

      // Update the model
      model->update();

    } catch (GRBException& e) {
      std::cerr << "Gurobi error: " << e.getMessage() << std::endl;
      throw std::runtime_error("Gurobi optimization failed");
    }
  }

  void initialize_model_three_reference_points() {
    try {
      const std::vector<int64_t>& r1 = refs[0];
      const std::vector<int64_t>& r2 = refs[1];
      const std::vector<int64_t>& r3 = refs[2];

      // ---------- Compute union / pairwise union reference points ----------
      std::vector<int64_t> r12(M), r13(M), r23(M), r123(M);
      for (int i = 0; i < M; ++i) {
        budget.check();
        r12[i] = std::max(r1[i], r2[i]);
        r13[i] = std::max(r1[i], r3[i]);
        r23[i] = std::max(r2[i], r3[i]);
        r123[i] = std::max({r1[i], r2[i], r3[i]});
      }

      // Build S for singles, pairs and triple
      GRBVar S1 = build_sum_of_logs(r1, "S_r1");
      GRBVar S2 = build_sum_of_logs(r2, "S_r2");
      GRBVar S3 = build_sum_of_logs(r3, "S_r3");
      GRBVar S12 = build_sum_of_logs(r12, "S_r12");
      GRBVar S13 = build_sum_of_logs(r13, "S_r13");
      GRBVar S23 = build_sum_of_logs(r23, "S_r23");
      GRBVar S123 = build_sum_of_logs(r123, "S_r123");

      // ---------- Compute LSEs for positives and negatives ----------

      // Inclusion-exclusion: union = (sum singles) - (sum pairs) + (triple)
      // We group positives = {S1, S2, S3, S123}; negatives = {S12, S13, S23}.
      std::vector<GRBVar> positives = {S1, S2, S3, S123};
      std::vector<GRBVar> negatives = {S12, S13, S23};

      GRBVar Lpos = log_sum_exp(positives, "pos_lse");
      GRBVar Lneg = log_sum_exp(negatives, "neg_lse");

      // Enforce Lpos > Lneg + EPS to guarantee 1 - exp(Lneg-Lpos) > 0
      model->addConstr(Lpos >= Lneg + EPS, "ensure_Lpos_gt_Lneg");

      // Compute B_minus_A = Lneg - Lpos  (<= -EPS)
      GRBVar B_minus_A = model->addVar(-GRB_INFINITY, 0.0, 0.0, GRB_CONTINUOUS, "B_minus_A");
      model->addConstr(B_minus_A == Lneg - Lpos, "def_B_minus_A");

      // exp(B_minus_A)
      GRBVar exp_B_minus_A = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "exp_B_minus_A");
      model->addGenConstrExp(B_minus_A, exp_B_minus_A, "exp_B_minus_A_constr");

      // one_minus = 1 - exp(B_minus_A)  (in (0,1])
      GRBVar one_minus = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "one_minus");
      model->addConstr(one_minus == 1.0 - exp_B_minus_A, "def_one_minus");

      // log_one_minus = log(one_minus)
      GRBVar log_one_minus = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "log_one_minus");
      model->addGenConstrLog(one_minus, log_one_minus, "log_one_minus_constr");

      // Final result: result = Lpos + log_one_minus  == log( exp(Lpos) - exp(Lneg) )
      GRBVar result = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "result");
      model->addConstr(result == Lpos + log_one_minus, "def_result");

      // Set objective
      model->setObjective(GRBLinExpr(result), GRB_MAXIMIZE);

      // Update model
      model->update();

    } catch (GRBException& e) {
      std::cerr << "Gurobi error: " << e.getMessage() << std::endl;
      throw std::runtime_error("Gurobi optimization failed");
    }
  }

  void initialize_model_general_reference_points() {
    try {
      int k = static_cast<int>(refs.size());

      // Inclusion–exclusion over all non-empty subsets of reference points
      std::vector<GRBVar> positives;
      std::vector<GRBVar> negatives;

      int total_masks = (1 << k);
      for (int mask = 1; mask < total_masks; ++mask) {
        budget.check();
        // Build union reference point for this subset
        std::vector<int64_t> rS(M, std::numeric_limits<int64_t>::min());
        for (int ref_idx = 0; ref_idx < k; ++ref_idx) {
          budget.check();
          if (mask & (1 << ref_idx)) {
            const std::vector<int64_t>& rj = refs[ref_idx];
            for (int i = 0; i < M; ++i) {
              budget.check();
              rS[i] = std::max(rS[i], rj[i]);
            }
          }
        }
        GRBVar S_mask = build_sum_of_logs(rS, "S_mask_" + std::to_string(mask));
        int bits = std::popcount(static_cast<unsigned int>(mask));
        if (bits % 2 == 1) {
          positives.push_back(S_mask);
        } else {
          negatives.push_back(S_mask);
        }
      }

      GRBVar Lpos = log_sum_exp(positives, "pos_lse");
      GRBVar Lneg = log_sum_exp(negatives, "neg_lse");

      // Enforce Lpos > Lneg + EPS to guarantee 1 - exp(Lneg-Lpos) > 0
      model->addConstr(Lpos >= Lneg + EPS, "ensure_Lpos_gt_Lneg");

      // Compute B_minus_A = Lneg - Lpos  (<= -EPS)
      GRBVar B_minus_A = model->addVar(-GRB_INFINITY, 0.0, 0.0, GRB_CONTINUOUS, "B_minus_A");
      model->addConstr(B_minus_A == Lneg - Lpos, "def_B_minus_A");

      // exp(B_minus_A)
      GRBVar exp_B_minus_A = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "exp_B_minus_A");
      model->addGenConstrExp(B_minus_A, exp_B_minus_A, "exp_B_minus_A_constr");

      // one_minus = 1 - exp(B_minus_A)  (in (0,1])
      GRBVar one_minus = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "one_minus");
      model->addConstr(one_minus == 1.0 - exp_B_minus_A, "def_one_minus");

      // log_one_minus = log(one_minus)
      GRBVar log_one_minus = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "log_one_minus");
      model->addGenConstrLog(one_minus, log_one_minus, "log_one_minus_constr");

      // Final result: result = Lpos + log_one_minus  == log( exp(Lpos) - exp(Lneg) )
      GRBVar result = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "result");
      model->addConstr(result == Lpos + log_one_minus, "def_result");

      // Set objective
      model->setObjective(GRBLinExpr(result), GRB_MAXIMIZE);

      // Update model
      model->update();

    } catch (GRBException& e) {
      std::cerr << "Gurobi error during model construction: " << e.getMessage() << std::endl;
      throw;
    }
  }

  // Get the model's number of variables and constraints
  std::tuple<int64_t, int64_t, std::vector<int64_t>> statistics() const {
    if (!model) return {0, 0, {0, 0}};
    int64_t num_vars = model->get(GRB_IntAttr_NumVars);
    int64_t num_lin_constrs = model->get(GRB_IntAttr_NumConstrs);
    int64_t num_qua_constrs = model->get(GRB_IntAttr_NumQConstrs);
    int64_t num_constrs = num_lin_constrs + num_qua_constrs;
    std::vector<int64_t> constr_types = {num_lin_constrs, num_qua_constrs};
    return std::make_tuple(num_vars, num_constrs, constr_types);
  }

  MOKPSolution solve() {
    if (initial_result.status != SolverStatus::NotSolved) return initial_result;
    MOKPSolution result;
    try {
      budget.check();
      model->set(GRB_DoubleParam_TimeLimit, budget.left());
      model->optimize();
      result.backend_status = model->get(GRB_IntAttr_Status);
      switch (result.backend_status) {
        case GRB_OPTIMAL: result.status = SolverStatus::Optimal; break;
        case GRB_INFEASIBLE: result.status = SolverStatus::Infeasible; break;
        case GRB_TIME_LIMIT: result.status = SolverStatus::TimeLimit; break;
        case GRB_INTERRUPTED: result.status = SolverStatus::Interrupted; break;
        case GRB_NUMERIC: result.status = SolverStatus::NumericalFailure; break;
        case GRB_UNBOUNDED: result.status = SolverStatus::Unbounded; break;
        default: result.status = SolverStatus::Other; break;
      }
      result.nodes = static_cast<int64_t>(model->get(GRB_DoubleAttr_NodeCount));
      if (model->get(GRB_IntAttr_SolCount) == 0) return result;
      result.items.assign(N, 0);
      result.values.assign(M, 0);
      for (int i = 0; i < N; ++i) {
        result.items[items[i].idx] = x[i].get(GRB_DoubleAttr_X) > 0.5;
        if (!result.items[items[i].idx]) continue;
        result.weight += items[i].weight;
        for (int d = 0; d < M; ++d) result.values[d] += items[i].values[d];
      }
      result.hypervolume = static_cast<hv_t>(std::round(std::exp(model->get(GRB_DoubleAttr_ObjVal))));
    } catch (GRBException& e) {
      result.status = SolverStatus::Error;
      result.message = e.getMessage();
    } catch (const TimeLimitReached& e) {
      result.status = SolverStatus::TimeLimit;
      result.message = e.what();
    } catch (const std::exception& e) {
      result.status = SolverStatus::NumericalFailure;
      result.message = e.what();
    }
    return result;
  }


};

} // namespace moco::ghps

#endif  // LOG_HPP
