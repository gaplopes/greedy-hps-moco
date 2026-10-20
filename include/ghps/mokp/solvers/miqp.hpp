#ifndef MIQP_HPP
#define MIQP_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ghps/mokp/instance.hpp>
#include <ghps/mokp/item.hpp>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "gurobi_c++.h"
#include <ghps/core/timer.hpp>

namespace moco::ghps {


/**
 * Hypervolume scalarization for one or more reference points.
 * Inclusion-exclusion combines products of integer gains f_i(x) - r_i.
 * Feasibility requires f_i(x) >= max_r(r_i) + 1 in every objective.
 * Normalized products use auxiliary variables and bilinear constraints with
 * NonConvex=2; the model objective is linear in those product variables.
 */

struct MOKPMIQPOptions {
  double time_limit = 600.0;
  std::vector<double> norm_vec = {};  // empty = automatic normalization for MOKP, ones for raw items
  int threads = 1;
  Timer deadline = {};
};

struct MOKPMIQP {
 private:
  Timer budget;
  MOKPSolution initial_result;
  std::unique_ptr<GRBEnv> env;
  std::unique_ptr<GRBModel> model;
  std::vector<Item> items;
  int32_t N;  // items
  int32_t M;  // objectives
  int64_t W;  // capacity

  std::vector<GRBVar> x;           // decision vars (binary)
  GRBConstr capacity_constr;       // Knapsack constraint
  std::vector<GRBVar> obj_values;  // f^i(x)

  // References
  std::vector<std::vector<int64_t>> refs;
  std::vector<int64_t> max_ref;

  // Normalization
  std::vector<double> norm_vec;

  // Model parameters
  double time_limit = 600.0;

  // Utility: create a var equal to a linear expression
  GRBVar var_from_expr(const GRBLinExpr& e, const std::string& name) {
    GRBVar v = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
    model->addConstr(v == e, name + "_def");
    return v;
  }

  // Utility: build chained product of linear expressions terms[0]*...*terms[k-1]
  // using k-1 bilinear equalities, returning the product variable
  GRBVar chained_product_of_terms(const std::vector<GRBLinExpr>& terms, const std::string& baseName) {
    if (terms.empty()) {
      // Create a constant-1 variable if called with empty vector
      GRBVar one = model->addVar(1.0, 1.0, 0.0, GRB_CONTINUOUS, baseName + "_one");
      return one;
    }
    // p0 = terms[0]
    GRBVar p_prev = var_from_expr(terms[0], baseName + "_p0");
    for (size_t k = 1; k < terms.size(); ++k) {
      budget.check();
      GRBVar p_k = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                                 baseName + "_p" + std::to_string(k));
      model->addQConstr(p_k == p_prev * terms[k], baseName + "_q" + std::to_string(k));
      p_prev = p_k;
    }
    return p_prev;  // final product
  }

  // Utility: build product of linear expressions using pairwise bilinear equalities.
  // For 3 terms (a,b,c): computes ab = a*b, then P = ab*c.
  // For 4 terms (a,b,c,d): computes ab = a*b, cd = c*d, then P = ab*cd.
  // For more terms, reduces by pairs level-by-level until one product remains.
  GRBVar product_of_terms(const std::vector<GRBLinExpr>& terms, const std::string& baseName) {
    if (terms.empty()) {
      // Create a constant-1 variable if called with empty vector
      GRBVar one = model->addVar(1.0, 1.0, 0.0, GRB_CONTINUOUS, baseName + "_one");
      return one;
    }

    // First, tie each linear term to its own variable so we can multiply variables pairwise
    std::vector<GRBVar> partials;
    partials.reserve(terms.size());
    for (size_t i = 0; i < terms.size(); ++i) {
      budget.check();
      GRBVar vi = var_from_expr(terms[i], baseName + "_t" + std::to_string(i));
      partials.push_back(vi);
    }

    // Reduce by pairwise multiplication: (a*b), (c*d), carry last if odd; repeat till one remains
    int level = 0;
    while (partials.size() > 1) {
      budget.check();
      std::vector<GRBVar> next;
      next.reserve((partials.size() + 1) / 2);
      for (size_t i = 0; i < partials.size(); i += 2) {
        budget.check();
        if (i + 1 < partials.size()) {
          std::string vname = baseName + "_L" + std::to_string(level) + "_p" + std::to_string(i / 2);
          std::string qname = baseName + "_L" + std::to_string(level) + "_q" + std::to_string(i / 2);
          GRBVar pij = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, vname);
          model->addQConstr(pij == partials[i] * partials[i + 1], qname);
          next.push_back(pij);
        } else {
          // carry over the last one
          next.push_back(partials[i]);
        }
      }
      partials.swap(next);
      ++level;
    }

    return partials.front();
  }

 public:
  MOKPMIQP(const MOKP& instance,
           const std::vector<std::vector<int64_t>>& refs,
           MOKPMIQPOptions options)
      : MOKPMIQP(instance.is_maximization ? instance.items : throw std::invalid_argument("MIQP supports maximization only"),
                 instance.N, instance.M, instance.W, refs, std::move(options),
                 &instance.dantzig_upper_bound()) {}

  /// Computes norm[i] = ub[i] / (abs_tol * product(ub))^(1/M).
  static std::vector<double> precision_optimal_norm(const std::vector<double>& upper_bound,
                                                  double abs_tol = 1e-6,
                                                  const Timer& deadline = Timer()) {
    deadline.check();
    if (std::any_of(upper_bound.begin(), upper_bound.end(), [](double value) { return value <= 0; }))
      return std::vector<double>(upper_bound.size(), 1.0);
    double log_product = std::log(abs_tol);
    for (double value : upper_bound) log_product += std::log(value);
    const double scale = std::exp(log_product / upper_bound.size());
    auto result = upper_bound;
    for (auto& value : result) value /= scale;
    deadline.check();
    return result;
  }

  MOKPMIQP(const std::vector<Item>& items, int32_t N, int32_t M, int64_t W,
           const std::vector<std::vector<int64_t>>& refs, MOKPMIQPOptions options)
      : MOKPMIQP(items, N, M, W, refs, std::move(options), nullptr) {}

 private:
  MOKPMIQP(const std::vector<Item>& items, int32_t N, int32_t M, int64_t W,
           const std::vector<std::vector<int64_t>>& refs, MOKPMIQPOptions options,
           const std::vector<double>* upper_bound)
      : budget(options.time_limit, options.deadline), items(items), N(N), M(M), W(W), refs(refs),
        time_limit(options.time_limit) {
    try {
      budget.check();

      norm_vec = std::move(options.norm_vec);
      if (norm_vec.empty()) norm_vec = upper_bound
          ? precision_optimal_norm(*upper_bound, 1e-6, budget) : std::vector<double>(M, 1.0);
      max_ref.assign(M, std::numeric_limits<int64_t>::min());
      for (const auto& r : refs)
        for (int i = 0; i < M; ++i) max_ref[i] = std::max(max_ref[i], r[i]);

      // Environment & model
      env = std::make_unique<GRBEnv>(true);
      env->set(GRB_IntParam_OutputFlag, 0);
      env->set(GRB_IntParam_Threads, options.threads);
      env->start();

      model = std::make_unique<GRBModel>(*env);
      model->set(GRB_DoubleParam_TimeLimit, time_limit);
      model->set(GRB_DoubleParam_MIPGap, 1e-6);
      model->set(GRB_IntParam_NumericFocus, 3);
      model->set(GRB_IntParam_DualReductions, 0);  // Resolve infeasible status without an input precheck.
      model->set(GRB_IntParam_NonConvex, 2);  // Enable bilinear equalities
      model->set(GRB_DoubleParam_FeasibilityTol, 1e-9);
      model->set(GRB_DoubleParam_IntFeasTol, 1e-9);



      // Decision vars: binary selection
      x.resize(N);
      for (int j = 0; j < N; ++j) {
        budget.check();
        x[j] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + std::to_string(j));
      }

      // Capacity
      GRBLinExpr capacity = 0;
      for (int j = 0; j < N; ++j) capacity += items[j].weight * x[j];
      capacity_constr = model->addConstr(capacity <= W, "knapsack_capacity");

      // Build f^i(x) values
      std::vector<GRBLinExpr> f_expr(M);
      obj_values.resize(M);
      for (int i = 0; i < M; ++i) {
        budget.check();
        f_expr[i] = 0;
        for (int j = 0; j < N; ++j) f_expr[i] += items[j].values[i] * x[j];
        obj_values[i] = var_from_expr(f_expr[i], "f_" + std::to_string(i));
      }

      // Integer strict dominance: f^i(x) >= max_ref^i + 1
      for (int i = 0; i < M; ++i) {
        budget.check();
        model->addConstr(obj_values[i] >= static_cast<double>(hv_t(max_ref[i]) + 1), "f_ge_max_ref_" + std::to_string(i));
      }

      // Feasibility constraints
      if (refs.size() == 1) {
        GRBVar P1;
        const auto& r = refs[0];

        // Build product P1 = prod_i ((f_i - r_i) / norm_vec_i)
        std::vector<GRBLinExpr> terms;
        terms.reserve(M);
        for (int i = 0; i < M; ++i) {
          budget.check();
          double inv = 1.0 / norm_vec[i];
          terms.emplace_back(obj_values[i] * inv - static_cast<double>(r[i]) * inv);
        }
        P1 = product_of_terms(terms, "P1");

        model->setObjective(GRBLinExpr(P1), GRB_MAXIMIZE);

      } else if (refs.size() == 2) {
        GRBVar P1, P2, PU;
        const auto& r1 = refs[0];
        const auto& r2 = refs[1];

        // P1 = prod_i (f_i - r1_i)
        {
          std::vector<GRBLinExpr> terms;
          terms.reserve(M);
          for (int i = 0; i < M; ++i) {
            budget.check();
            double inv = 1.0 / norm_vec[i];
            terms.emplace_back(obj_values[i] * inv - static_cast<double>(r1[i]) * inv);
          }
          P1 = product_of_terms(terms, "P1");
        }

        // P2 = prod_i (f_i - r2_i)
        {
          std::vector<GRBLinExpr> terms;
          terms.reserve(M);
          for (int i = 0; i < M; ++i) {
            budget.check();
            double inv = 1.0 / norm_vec[i];
            terms.emplace_back(obj_values[i] * inv - static_cast<double>(r2[i]) * inv);
          }
          P2 = product_of_terms(terms, "P2");
        }

        // PU = prod_i (f_i - max(r1_i, r2_i))
        {
          std::vector<GRBLinExpr> terms;
          terms.reserve(M);
          for (int i = 0; i < M; ++i) {
            budget.check();
            double inv = 1.0 / norm_vec[i];
            terms.emplace_back(obj_values[i] * inv - static_cast<double>(std::max(r1[i], r2[i])) * inv);
          }
          PU = product_of_terms(terms, "PU");
        }

        // Objective: P1 + P2 - PU
        model->setObjective(GRBLinExpr(P1) + GRBLinExpr(P2) - GRBLinExpr(PU), GRB_MAXIMIZE);
      } else if (refs.size() == 3) {
        GRBVar P1, P2, P3, PU12, PU13, PU23, PU123;
        // Inclusion-Exclusion for k=3 reference points
        const auto& r1 = refs[0];
        const auto& r2 = refs[1];
        const auto& r3 = refs[2];

        // Individual products: P1, P2, P3
        {
          std::vector<GRBLinExpr> t1, t2, t3;
          t1.reserve(M);
          t2.reserve(M);
          t3.reserve(M);
          for (int i = 0; i < M; ++i) {
            budget.check();
            double inv = 1.0 / norm_vec[i];
            t1.emplace_back(obj_values[i] * inv - static_cast<double>(r1[i]) * inv);
            t2.emplace_back(obj_values[i] * inv - static_cast<double>(r2[i]) * inv);
            t3.emplace_back(obj_values[i] * inv - static_cast<double>(r3[i]) * inv);
          }
          P1 = product_of_terms(t1, "P1");
          P2 = product_of_terms(t2, "P2");
          P3 = product_of_terms(t3, "P3");
        }

        // Pairwise max products: PU12, PU13, PU23
        {
          std::vector<GRBLinExpr> u12, u13, u23;
          u12.reserve(M);
          u13.reserve(M);
          u23.reserve(M);
          for (int i = 0; i < M; ++i) {
            budget.check();
            double inv = 1.0 / norm_vec[i];
            u12.emplace_back(obj_values[i] * inv - static_cast<double>(std::max(r1[i], r2[i])) * inv);
            u13.emplace_back(obj_values[i] * inv - static_cast<double>(std::max(r1[i], r3[i])) * inv);
            u23.emplace_back(obj_values[i] * inv - static_cast<double>(std::max(r2[i], r3[i])) * inv);
          }
          PU12 = product_of_terms(u12, "PU12");
          PU13 = product_of_terms(u13, "PU13");
          PU23 = product_of_terms(u23, "PU23");
        }

        // PU123 = prod_i (f_i - max(r1_i, r2_i, r3_i))
        {
          std::vector<GRBLinExpr> u123;
          u123.reserve(M);
          for (int i = 0; i < M; ++i) {
            budget.check();
            double inv = 1.0 / norm_vec[i];
            int64_t m123 = std::max(r1[i], std::max(r2[i], r3[i]));
            u123.emplace_back(obj_values[i] * inv - static_cast<double>(m123) * inv);
          }
          PU123 = product_of_terms(u123, "PU123");
        }

        // (P1 + P2 + P3) - (PU12 + PU13 + PU23) + PU123
        GRBLinExpr objExpr = GRBLinExpr(P1) + GRBLinExpr(P2) + GRBLinExpr(P3) - GRBLinExpr(PU12) - GRBLinExpr(PU13) - GRBLinExpr(PU23) + GRBLinExpr(PU123);
        model->setObjective(objExpr, GRB_MAXIMIZE);
      } else {
        // General inclusion–exclusion for arbitrary k >= 4 reference points
        const int k = static_cast<int>(refs.size());

        // Build inclusion–exclusion objective:
        // Sum over all non-empty subsets S of {1..k} with sign (-1)^{|S|+1}
        GRBLinExpr objExpr = 0;
        const int totalMasks = (1 << k);
        for (int mask = 1; mask < totalMasks; ++mask) {
          budget.check();
          // Compute per-objective max reference over subset S
          std::vector<GRBLinExpr> terms;
          terms.reserve(M);
          for (int i = 0; i < M; ++i) {
            budget.check();
            int64_t max_ref = std::numeric_limits<int64_t>::min();
            for (int j = 0; j < k; ++j) {
              budget.check();
              if (mask & (1 << j)) max_ref = std::max<int64_t>(max_ref, refs[j][i]);
            }
            double inv = 1.0 / norm_vec[i];
            terms.emplace_back(obj_values[i] * inv - static_cast<double>(max_ref) * inv);
          }

          // Product for this subset
          GRBVar P_S = product_of_terms(terms, "PS_" + std::to_string(mask));

          // Inclusion–exclusion sign
          int bits = std::popcount(static_cast<unsigned int>(mask));
          int sign = (bits % 2 == 1) ? +1 : -1;
          objExpr += sign * GRBLinExpr(P_S);
        }

        model->setObjective(objExpr, GRB_MAXIMIZE);
      }

      model->update();
    } catch (GRBException& e) {
      std::cerr << "Gurobi error during model construction: " << e.getMessage() << std::endl;
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

 public:
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
      double hypervolume = model->get(GRB_DoubleAttr_ObjVal);
      for (double norm : norm_vec) hypervolume *= norm;
      result.hypervolume = static_cast<hv_t>(std::round(hypervolume));
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

#endif  // MIQP_HPP
