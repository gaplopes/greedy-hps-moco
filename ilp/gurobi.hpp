#ifndef GUROBI_HPP
#define GUROBI_HPP

#include <structs.hpp>
#include "gurobi_c++.h"

// --------------------------------------------

inline void solve_3D_3R_QUAD(GRBModel& model, const ProblemVariables& problem, std::vector<RefPoint>& refPoints, int64_t& objective_value, std::string& solution_items, std::vector<int64_t>& solution_point) {
  // Create variables
  const int32_t N = problem.numElements;
  const int32_t DIM = problem.numObjectives;

  // Get reference points
  std::vector<int64_t> refPointOne = refPoints[0].coordinates;
  std::vector<int64_t> refPointTwo = refPoints[1].coordinates;
  std::vector<int64_t> refPointThree = refPoints[2].coordinates;
  std::vector<int64_t> refPointInterOneTwo = std::vector<int64_t>(DIM);
  for (int i = 0; i < DIM; i++) {
    refPointInterOneTwo[i] = std::max(refPointOne[i], refPointTwo[i]);
  }
  std::vector<int64_t> refPointInterOneThree = std::vector<int64_t>(DIM);
  for (int i = 0; i < DIM; i++) {
    refPointInterOneThree[i] = std::max(refPointOne[i], refPointThree[i]);
  }
  std::vector<int64_t> refPointInterTwoThree = std::vector<int64_t>(DIM);
  for (int i = 0; i < DIM; i++) {
    refPointInterTwoThree[i] = std::max(refPointTwo[i], refPointThree[i]);
  }
  std::vector<int64_t> refPointInterOneTwoThree = std::vector<int64_t>(DIM);
  for (int i = 0; i < DIM; i++) {
    refPointInterOneTwoThree[i] = std::max(refPointOne[i], std::max(refPointTwo[i], refPointThree[i]));
  }

  // Create GUROBI variables
  GRBVar y_ij[N][N], x_k[N];
  for (int32_t i = 0; i < N; i++) {
    for (int32_t j = 0; j < N; j++) {
      y_ij[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
    x_k[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
  }

  model.update();

  // Set objective
  GRBQuadExpr objExpr = 0;

  for (int32_t i = 0; i < N; i++) {
    for (int32_t j = 0; j < N; j++) {
      for (int32_t k = 0; k < N; k++) {
        objExpr += problem.values[i][0] * problem.values[j][1] * problem.values[k][2] * y_ij[i][j] * x_k[k];
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      auto compute_value = [&problem, i, j](std::vector<int64_t> refPoint) {
        int64_t value = (problem.values[i][0] * problem.values[j][2] * refPoint[1]) +
                        (problem.values[i][1] * problem.values[j][2] * refPoint[0]) +
                        (problem.values[i][0] * problem.values[j][1] * refPoint[2]);
        return value;
      };
      int64_t value = compute_value(refPointOne) + compute_value(refPointTwo) + compute_value(refPointThree) -
                      compute_value(refPointInterOneTwo) - compute_value(refPointInterOneThree) -
                      compute_value(refPointInterTwoThree) + compute_value(refPointInterOneTwoThree);
      objExpr -= value * y_ij[i][j];
    }
  }

  for (int i = 0; i < N; i++) {
    auto compute_value = [&problem, i](std::vector<int64_t> refPoint) {
      int64_t value = (problem.values[i][0] * refPoint[1] * refPoint[2]) +
                      (problem.values[i][1] * refPoint[0] * refPoint[2]) +
                      (problem.values[i][2] * refPoint[0] * refPoint[1]);
      return value;
    };
    int64_t value = compute_value(refPointOne) + compute_value(refPointTwo) + compute_value(refPointThree) -
                    compute_value(refPointInterOneTwo) - compute_value(refPointInterOneThree) -
                    compute_value(refPointInterTwoThree) + compute_value(refPointInterOneTwoThree);
    objExpr += value * y_ij[i][i];
  }

  int64_t constantTermRefOne = (refPointOne[0] * refPointOne[1] * refPointOne[2]);
  int64_t constantTermRefTwo = (refPointTwo[0] * refPointTwo[1] * refPointTwo[2]);
  int64_t constantTermRefThree = (refPointThree[0] * refPointThree[1] * refPointThree[2]);
  int64_t constantTermRefInterOneTwo = (refPointInterOneTwo[0] * refPointInterOneTwo[1] * refPointInterOneTwo[2]);
  int64_t constantTermRefInterOneThree = (refPointInterOneThree[0] * refPointInterOneThree[1] * refPointInterOneThree[2]);
  int64_t constantTermRefInterTwoThree = (refPointInterTwoThree[0] * refPointInterTwoThree[1] * refPointInterTwoThree[2]);
  int64_t constantTermRefInterOneTwoThree = (refPointInterOneTwoThree[0] * refPointInterOneTwoThree[1] * refPointInterOneTwoThree[2]);
  int64_t constantTerm = constantTermRefOne + constantTermRefTwo + constantTermRefThree - constantTermRefInterOneTwo -
                         constantTermRefInterOneThree - constantTermRefInterTwoThree + constantTermRefInterOneTwoThree;
  objExpr -= constantTerm;

  model.setObjective(objExpr, GRB_MAXIMIZE);

  // Add constraints

  GRBQuadExpr weightConstraint = 0;
  // Weight constraint: \sum_{i=1}^{n} w_{i} y_{ii}x_{i} <= W
  for (int i = 0; i < N; i++) {
    weightConstraint += problem.weights[i] * y_ij[i][i] * x_k[i];
  }
  model.addQConstr(weightConstraint <= problem.maxWeight, "Constraint_W");

  GRBQuadExpr refPointConstraint1 = 0, refPointConstraint2 = 0, refPointConstraint3 = 0;
  // Ref point constraints: \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{1}
  //                        \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{2}
  //                        \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{3}
  for (int i = 0; i < N; i++) {
    refPointConstraint1 += problem.values[i][0] * y_ij[i][i] * x_k[i];
    refPointConstraint2 += problem.values[i][1] * y_ij[i][i] * x_k[i];
    refPointConstraint3 += problem.values[i][2] * y_ij[i][i] * x_k[i];
  }
  model.addQConstr(refPointConstraint1 >= refPointInterOneTwoThree[0], "Constraint_R1");
  model.addQConstr(refPointConstraint2 >= refPointInterOneTwoThree[1], "Constraint_R2");
  model.addQConstr(refPointConstraint3 >= refPointInterOneTwoThree[2], "Constraint_R3");

  // Other constraints: y_{ii} x_{i} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{jj} x_{j} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{kk} x_{k} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{ij} == x_{i} x_{j}           \forall i,j \in N
  //                    y_{ii}x_{i} + y_{jj}x_{j} + y_{kk}x_{k} - 2 <= y_{ij} x_{k}     \forall i,j,k \in N
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        model.addQConstr(y_ij[i][i] * x_k[i] >= y_ij[i][j] * x_k[k], "Constraint_1");
        model.addQConstr(y_ij[j][j] * x_k[j] >= y_ij[i][j] * x_k[k], "Constraint_2");
        model.addQConstr(y_ij[k][k] * x_k[k] >= y_ij[i][j] * x_k[k], "Constraint_3");
        model.addQConstr(y_ij[i][i] * x_k[i] + y_ij[j][j] * x_k[j] + y_ij[k][k] * x_k[k] - 2 <= y_ij[i][j] * x_k[k], "Constraint_5");
      }
      model.addQConstr(y_ij[i][j] == x_k[i] * x_k[j], "Constraint_4");
    }
  }

  // Optimize model
  model.optimize();

  // Verify if time limit was reached or model is infeasible
  int32_t status = model.get(GRB_IntAttr_Status);
  if (status == GRB_TIME_LIMIT || status == GRB_INFEASIBLE) {
    return;
  }

  // Get objective value
  objective_value = model.get(GRB_DoubleAttr_ObjVal);

  // Get solution
  solution_items = "";
  solution_point = std::vector<int64_t>(DIM, 0);
  for (int32_t i = 0; i < N; i++) {
    if (x_k[i].get(GRB_DoubleAttr_X) > 0.5) {
      solution_items += "1";
      for (int32_t j = 0; j < DIM; j++) {
        solution_point[j] += problem.values[i][j];
      }
    } else {
      solution_items += "0";
    }
    solution_items += (i < N - 1) ? " " : "";
  }
}

inline void solve_3D_2R_QUAD(GRBModel& model, const ProblemVariables& problem, std::vector<RefPoint>& refPoints, int64_t& objective_value, std::string& solution_items, std::vector<int64_t>& solution_point) {
  // Create variables
  const int32_t N = problem.numElements;
  const int32_t DIM = problem.numObjectives;

  // Get reference points
  std::vector<int64_t> refPointOne = refPoints[0].coordinates;
  std::vector<int64_t> refPointTwo = refPoints[1].coordinates;
  std::vector<int64_t> refPointInter(DIM);
  for (int i = 0; i < DIM; i++) {
    refPointInter[i] = std::max(refPointOne[i], refPointTwo[i]);
  }

  // Create GUROBI variables
  GRBVar y_ij[N][N], x_k[N];
  for (int32_t i = 0; i < N; i++) {
    for (int32_t j = 0; j < N; j++) {
      y_ij[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
    x_k[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
  }

  model.update();

  // Set objective
  GRBQuadExpr objExpr = 0;

  for (int32_t i = 0; i < N; i++) {
    for (int32_t j = 0; j < N; j++) {
      for (int32_t k = 0; k < N; k++) {
        objExpr += problem.values[i][0] * problem.values[j][1] * problem.values[k][2] * y_ij[i][j] * x_k[k];
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      auto compute_value = [&problem, i, j](std::vector<int64_t> refPoint) {
        int64_t value = (problem.values[i][0] * problem.values[j][2] * refPoint[1]) +
                        (problem.values[i][1] * problem.values[j][2] * refPoint[0]) +
                        (problem.values[i][0] * problem.values[j][1] * refPoint[2]);
        return value;
      };
      objExpr -= (compute_value(refPointOne) + compute_value(refPointTwo) - compute_value(refPointInter)) * y_ij[i][j];
    }
  }

  for (int i = 0; i < N; i++) {
    auto compute_value = [&problem, i](std::vector<int64_t> refPoint) {
      int64_t value = (problem.values[i][0] * refPoint[1] * refPoint[2]) +
                      (problem.values[i][1] * refPoint[0] * refPoint[2]) +
                      (problem.values[i][2] * refPoint[0] * refPoint[1]);
      return value;
    };
    objExpr += (compute_value(refPointOne) + compute_value(refPointTwo) - compute_value(refPointInter)) * y_ij[i][i];
  }

  int64_t constantTermRefOne = (refPointOne[0] * refPointOne[1] * refPointOne[2]);
  int64_t constantTermRefTwo = (refPointTwo[0] * refPointTwo[1] * refPointTwo[2]);
  int64_t constantTermRefInter = (refPointInter[0] * refPointInter[1] * refPointInter[2]);
  int64_t constantTerm = constantTermRefOne + constantTermRefTwo - constantTermRefInter;
  objExpr -= constantTerm;

  model.setObjective(objExpr, GRB_MAXIMIZE);

  // Add constraints

  GRBQuadExpr weightConstraint = 0;
  // Weight constraint: \sum_{i=1}^{n} w_{i} y_{ii}x_{i} <= W
  for (int i = 0; i < N; i++) {
    weightConstraint += problem.weights[i] * y_ij[i][i] * x_k[i];
  }
  model.addQConstr(weightConstraint <= problem.maxWeight, "Constraint_W");

  GRBQuadExpr refPointConstraint1 = 0, refPointConstraint2 = 0, refPointConstraint3 = 0;
  // Ref point constraints: \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{1}
  //                        \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{2}
  //                        \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{3}
  for (int i = 0; i < N; i++) {
    refPointConstraint1 += problem.values[i][0] * y_ij[i][i] * x_k[i];
    refPointConstraint2 += problem.values[i][1] * y_ij[i][i] * x_k[i];
    refPointConstraint3 += problem.values[i][2] * y_ij[i][i] * x_k[i];
  }
  model.addQConstr(refPointConstraint1 >= refPointInter[0], "Constraint_R1");
  model.addQConstr(refPointConstraint2 >= refPointInter[1], "Constraint_R2");
  model.addQConstr(refPointConstraint3 >= refPointInter[2], "Constraint_R3");

  // Other constraints: y_{ii} x_{i} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{jj} x_{j} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{kk} x_{k} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{ij} == x_{i} x_{j}           \forall i,j \in N
  //                    y_{ii}x_{i} + y_{jj}x_{j} + y_{kk}x_{k} - 2 <= y_{ij} x_{k}     \forall i,j,k \in N
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        model.addQConstr(y_ij[i][i] * x_k[i] >= y_ij[i][j] * x_k[k], "Constraint_1");
        model.addQConstr(y_ij[j][j] * x_k[j] >= y_ij[i][j] * x_k[k], "Constraint_2");
        model.addQConstr(y_ij[k][k] * x_k[k] >= y_ij[i][j] * x_k[k], "Constraint_3");
        model.addQConstr(y_ij[i][i] * x_k[i] + y_ij[j][j] * x_k[j] + y_ij[k][k] * x_k[k] - 2 <= y_ij[i][j] * x_k[k], "Constraint_5");
      }
      model.addQConstr(y_ij[i][j] == x_k[i] * x_k[j], "Constraint_4");
    }
  }

  // Optimize model
  model.optimize();

  // Verify if time limit was reached or model is infeasible
  int32_t status = model.get(GRB_IntAttr_Status);
  if (status == GRB_TIME_LIMIT || status == GRB_INFEASIBLE) {
    return;
  }

  // Get objective value
  objective_value = model.get(GRB_DoubleAttr_ObjVal);

  // Get solution
  solution_items = "";
  solution_point = std::vector<int64_t>(DIM, 0);
  for (int32_t i = 0; i < N; i++) {
    if (x_k[i].get(GRB_DoubleAttr_X) > 0.5) {
      solution_items += "1";
      for (int32_t j = 0; j < DIM; j++) {
        solution_point[j] += problem.values[i][j];
      }
    } else {
      solution_items += "0";
    }
    solution_items += (i < N - 1) ? " " : "";
  }
}

inline void solve_3D_1R_QUAD(GRBModel& model, const ProblemVariables& problem, std::vector<RefPoint>& refPoints, int64_t& objective_value, std::string& solution_items, std::vector<int64_t>& solution_point) {
  // Create variables
  const int32_t N = problem.numElements;
  const int32_t DIM = problem.numObjectives;

  std::vector<int64_t> refPoint = refPoints[0].coordinates;

  GRBVar y_ij[N][N], x_k[N];
  for (int32_t i = 0; i < N; i++) {
    for (int32_t j = 0; j < N; j++) {
      y_ij[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
    x_k[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
  }

  model.update();

  // Set objective
  GRBQuadExpr objExpr = 0;

  for (int32_t i = 0; i < N; i++) {
    for (int32_t j = 0; j < N; j++) {
      for (int32_t k = 0; k < N; k++) {
        objExpr += problem.values[i][0] * problem.values[j][1] * problem.values[k][2] * y_ij[i][j] * x_k[k];
      }
    }
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      auto compute_value = [&problem, i, j](std::vector<int64_t> refPoint) {
        int64_t value = (problem.values[i][0] * problem.values[j][2] * refPoint[1]) +
                        (problem.values[i][1] * problem.values[j][2] * refPoint[0]) +
                        (problem.values[i][0] * problem.values[j][1] * refPoint[2]);
        return value;
      };
      objExpr -= compute_value(refPoint) * y_ij[i][j];
    }
  }

  for (int i = 0; i < N; i++) {
    auto compute_value = [&problem, i](std::vector<int64_t> refPoint) {
      int64_t value = (problem.values[i][0] * refPoint[1] * refPoint[2]) +
                      (problem.values[i][1] * refPoint[0] * refPoint[2]) +
                      (problem.values[i][2] * refPoint[0] * refPoint[1]);
      return value;
    };
    objExpr += compute_value(refPoint) * y_ij[i][i];
  }

  objExpr -= (refPoint[0] * refPoint[1] * refPoint[2]);

  model.setObjective(objExpr, GRB_MAXIMIZE);

  // Add constraints

  GRBQuadExpr weightConstraint = 0;
  // Weight constraint: \sum_{i=1}^{n} w_{i} y_{ii}x_{i} <= W
  for (int i = 0; i < N; i++) {
    weightConstraint += problem.weights[i] * y_ij[i][i] * x_k[i];
  }
  model.addQConstr(weightConstraint <= problem.maxWeight, "Constraint_W");

  GRBQuadExpr refPointConstraint1 = 0, refPointConstraint2 = 0, refPointConstraint3 = 0;
  // Ref point constraints: \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{1}
  //                        \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{2}
  //                        \sum_{i=1}^{n} p_{i}^{j} y_{ii}x_{i} >= r^{3}
  for (int i = 0; i < N; i++) {
    refPointConstraint1 += problem.values[i][0] * y_ij[i][i] * x_k[i];
    refPointConstraint2 += problem.values[i][1] * y_ij[i][i] * x_k[i];
    refPointConstraint3 += problem.values[i][2] * y_ij[i][i] * x_k[i];
  }
  model.addQConstr(refPointConstraint1 >= refPoint[0], "Constraint_R1");
  model.addQConstr(refPointConstraint2 >= refPoint[1], "Constraint_R2");
  model.addQConstr(refPointConstraint3 >= refPoint[2], "Constraint_R3");

  // Other constraints: y_{ii} x_{i} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{jj} x_{j} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{kk} x_{k} >= y_{ij} x_{k}    \forall i,j,k \in N
  //                    y_{ij} == x_{i} x_{j}           \forall i,j \in N
  //                    y_{ii}x_{i} + y_{jj}x_{j} + y_{kk}x_{k} - 2 <= y_{ij} x_{k}     \forall i,j,k \in N
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        model.addQConstr(y_ij[i][i] * x_k[i] >= y_ij[i][j] * x_k[k], "Constraint_1");
        model.addQConstr(y_ij[j][j] * x_k[j] >= y_ij[i][j] * x_k[k], "Constraint_2");
        model.addQConstr(y_ij[k][k] * x_k[k] >= y_ij[i][j] * x_k[k], "Constraint_3");
        model.addQConstr(y_ij[i][i] * x_k[i] + y_ij[j][j] * x_k[j] + y_ij[k][k] * x_k[k] - 2 <= y_ij[i][j] * x_k[k], "Constraint_5");
      }
      model.addQConstr(y_ij[i][j] == x_k[i] * x_k[j], "Constraint_4");
    }
  }

  // Optimize model
  model.optimize();

  // Verify if time limit was reached or model is infeasible
  int32_t status = model.get(GRB_IntAttr_Status);
  if (status == GRB_TIME_LIMIT || status == GRB_INFEASIBLE) {
    return;
  }

  // Get objective value
  objective_value = model.get(GRB_DoubleAttr_ObjVal);

  // Get solution
  solution_items = "";
  solution_point = std::vector<int64_t>(DIM, 0);
  for (int32_t i = 0; i < N; i++) {
    if (x_k[i].get(GRB_DoubleAttr_X) > 0.5) {
      solution_items += "1";
      for (int32_t j = 0; j < DIM; j++) {
        solution_point[j] += problem.values[i][j];
      }
    } else {
      solution_items += "0";
    }
    solution_items += (i < N - 1) ? " " : "";
  }
}

auto solveILP(const ProblemVariables& problem, const std::vector<RefPoint>& refPoints, const int32_t DIM, const int32_t N_REF, const int32_t MAX_TIME, bool& finished) {
  int64_t objective = -1;
  std::string solution = "";
  std::vector<int64_t> solution_point(DIM, 0);

  try {
    // Create environment
    GRBEnv env = GRBEnv(true);
    // env.set("LogFile", "gurobi.log");              // Set log file
    env.set(GRB_IntParam_OutputFlag, 0);           // Disable output
    env.set(GRB_DoubleParam_TimeLimit, MAX_TIME);  // Set time limit in seconds
    env.set(GRB_DoubleParam_MIPGap, 1e-9);         // Set relative MIP optimality gap
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Set objective to maximize
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    model.set(GRB_IntParam_Method, 2);        // Barrier method
    model.set(GRB_IntParam_NumericFocus, 1);  // Numerical focus on optimality

    // --------------------------------------------

    switch (DIM) {
      case 3:
        switch (N_REF) {
          case 1:
            solve_3D_1R_QUAD(model, problem, refPoints, objective, solution, solution_point);
            break;
          case 2:
            solve_3D_2R_QUAD(model, problem, refPoints, objective, solution, solution_point);
            break;
          case 3:
            solve_3D_3R_QUAD(model, problem, refPoints, objective, solution, solution_point);
            break;
          default:
            break;
        }
        break;
      default:
        break;
    }

    // --------------------------------------------

    int32_t status = model.get(GRB_IntAttr_Status);
    if (status == GRB_TIME_LIMIT || status == GRB_INFEASIBLE) {
      finished = false;
    }

  } catch (GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch (...) {
    std::cout << "Exception during optimization" << std::endl;
  }

  // std::cout << "Objective: " << objective << std::endl;
  // std::cout << "Solution: " << solution << std::endl;
  // std::cout << "Solution point: ";
  // for (auto value : solution_point) {
  //   std::cout << value << " ";
  // }
  // std::cout << std::endl;

  return std::make_tuple(objective, solution, solution_point);
}

#endif  // GUROBI_HPP
