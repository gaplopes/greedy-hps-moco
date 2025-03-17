# greedy-hps-moco

A Greedy Hypervolume Polychotomic Scheme for Multiobjective Combinatorial Optimization (Greedy-HPS-MOCO).

## Introduction

Greedy-HPS-MOCO is an approach to solve Multiobjective Combinatorial Optimization (MOCO) problems.
It introduces a generic greedy method to compute representations of the nondominated set that approximately maximizes the dominated hypervolume.

In particular, the Greedy-HPS-MOCO is designed to solve the Multiobjective Knapsack Problem (MOKP), maximizing the hypervolume dominated by a 
subset of nondominated solutions of size J with respect to K reference points.

There are three solver variants available:

1. **Greedy**: The main solver that uses a greedy algorithm to find a representation set. In this case, the non-dominated set is required to be known.
2. **ILP**: An alternative solver that uses a mixed-integer programming formulation for hypervolume scalarization.
3. **BB**: A specialized branch-and-bound solver for m-objective knapsack problems.

In both the ILP and BB solvers, the non-dominated set is not required to be known, and the algorithm will compute it during the solution process.
Moreover, for both solvers, it is necessary to provide the number of reference points K to be used in the hypervolume scalarization.

## Key Features

- Hypervolume Scalarization: Iteratively builds a representation by solving a sequence of hypervolume scalarized problems
- Flexible Reference Points: Supports k reference points as a configurable parameter
- Mixed-Integer Formulation: Includes a mixed-integer programming formulation of the hypervolume scalarization function
- Specialized Branch-and-Bound: Implements a combinatorial branch-and-bound algorithm specifically tailored for m-objective knapsack problems

## Problem Context

Multiobjective optimization problems involve optimizing several conflicting objectives simultaneously.
The solution to such problems is typically not a single point but a set of Pareto-optimal solutions representing different trade-offs among the objectives.

This implementation focuses on:

1. Finding a representative subset of nondominated solutions
2. Maximizing the hypervolume dominated by this representation
3. Balancing computational efficiency with solution quality

## Requirements
- C++17 compatible compiler
- CMake 3.10 or higher
- Gurobi 9.0 or higher (optional, for mixed-integer formulation)
- CLI11 (included as a submodule)
- OpenMP (optional, for parallelization)

## Building

To build the Greedy-HPS-MOCO solver, follow these steps:

```bash
# Clone the repository
git clone
cd greedy-hps-moco
mkdir build
cd build
cmake ..
cmake --build .
```

## Usage

Basic usage of the Greedy-HPS-MOCO solver is as follows:

```bash
# Solve a 3-objective knapsack problem with 10 items, using 2 reference points
./greedy-hps-moco --input-file={instance_path} --J=10 --K=2 --bb=true --detailed-output=true
```

Key parameters include:

- `--input-file`: Path to the input file containing the instance data
- `--J`: Size of the representation set to find
- `--K`: Number of reference points to consider
- `--time-limit`: Maximum time limit for the solver in seconds (default: 3600)
- `--n-threads`: Number of threads to use for parallelization (default: 1)
- `--ilp`: Use the mixed-integer formulation for hypervolume scalarization (requires Gurobi)
- `--bb`: Use the specialized branch-and-bound algorithm for m-objective knapsack problems (default)
- `--greedy`: Use the greedy algorithm (only works if the non-dominated set is not empty)
- `--detailed-output`: Output detailed statistics and solution coordinates

## File Formats

The input file format expected for the MOKP instances is as follows:

```
N M # Number of items, number of objectives
W   # Knapsack capacity
w1 v11 v12 ... v1M # Item 1 weight, objective values
w2 v21 v22 ... v2M # Item 2 weight, objective values
...
wN vN1 vN2 ... vNM # Item N weight, objective values
```

The instances considered for the experiments are available in a repository at [https://github.com/gaplopes/mobkp-instances].

## Output Information

Greedy-HPS-MOCO provides comprehensive statistics about the solution process and quality of results.
The output includes performance metrics, solution coordinates, and hypervolume calculations.

### Solution Details

When using the detailed output mode, the solver provides:

- **Solution Coordinates**: The objective values for each solution in the representation set
- **Performance Statistics**: Iteration-by-iteration statistics in the format (k, subproblems, hypervolume)
- **Elapsed Time**: Total time taken to find the representation set

### Algorithm Performance Metrics

The output includes detailed performance indicators:

- **Iterations**: Total number of iterations performed by the algorithm
- **Subproblems**: Number of subproblems solved during the branch-and-bound process
- **Cache Hits**: Number of cache hits during the solution process
- **Time Breakdown**: 
  - Subgraphs computation time
  - Subproblems solving time
  - Total elapsed time

### Quality Metrics

The solver calculates various hypervolume-based quality metrics:

- **HV(Nset)**: Hypervolume dominated by the complete nondominated set (if available)
- **HV(Solutions)**: Hypervolume dominated by the computed representation set
- **HV Ratio**: The ratio between the representation set's hypervolume and the complete set's hypervolume
- **Nadir-based Metrics**: Alternative hypervolume calculations using the nadir point of the nondominated set
  - HV(Nadir Nset): Hypervolume of the nondominated set using the nadir point as reference
  - HV(Nadir Solutions): Hypervolume of the representation using the nadir point
  - HV Ratio (Nadir): The corresponding ratio with nadir-based calculations

### Set Comparison

- **|Nset|**: Size of the complete nondominated set (if available)
- **|Solutions|**: Size of the computed representation set
- **|Matching|**: Number of solutions that are common between the representation and the complete set

### Sample Output

```
Statistics:
(1,125,0.03) (2,250,1.02) (1,375,2.01) (4,200,0.01)
Iterations: 4
Subproblems: 128
Cache hits: 32
Subgraphs time: 0.125s
Subproblems time: 0.375s
Solutions:
(100,200,300) (400,300,200) (350,250,150) (250,350,150)
Elapsed time: 0.5s
|Nset|: 20
|Solutions|: 4
|Matching|: 2
HV(Nset): 5000000
HV(Solutions): 4500000
HV ratio: 0.9
HV(Nadir Nset): 4800000
HV(Nadir Solutions): 4300000
HV ratio (Nadir): 0.896
```

### Output Formats

The solver supports two output formats:
- **Detailed**: Human-readable format with labels and newlines (use `--detailed-output` flag)
- **Compact**: Space-separated values suitable for parsing by analysis scripts (default)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
