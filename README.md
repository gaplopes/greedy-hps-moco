# greedy-hps-moco

A Greedy Hypervolume Polychotomic Scheme for Multiobjective Combinatorial
Optimization (Greedy-HPS-MOCO), implemented as a C++20 header-only library for
multiobjective knapsack problems.

## Behavior

Knapsack objectives are maximized. GHPS builds a representation containing **up to
J solutions with positive hypervolume gain**. It stops when no positive gain
remains or the computation terminates. It does not enumerate the complete Pareto
front or certify a globally optimal size-J representation.

`K` is the maximum number of connected local reference bounds used in one
scalarization. The representation's total hypervolume is measured relative to the
instance's all-zero reference point.

The solver CLI selects its backend with `--mode`:

| Mode | Behavior | Requires a supplied nondominated set | Requires a runtime Gurobi license |
| --- | --- | --- | --- |
| `0` | Branch-and-bound scalarizations | No | No |
| `1` | Product model with bilinear constraints, named MIQP | No | Yes |
| `2` | Log-exp scalarizations | No | Yes |
| `3` | Log-exp for subsets of size at most two, otherwise MIQP | No | Yes |
| `4` | Greedy selection from the supplied front; default | Yes | No |

A complete supplied front is needed for trustworthy comparisons against
simulation. Parsing validates its shape and coordinate domain, not its
completeness or feasibility.

## Requirements

- A C++20 compiler supporting `__int128`.
- CMake 3.18 or newer.
- Gurobi SDK 12.0 or newer. Runtime checks have been run with Gurobi 12.0.3.
- OpenMP for parallel HPS subproblems.
- mooutils, found as a package or fetched at a pinned revision.
- CLI11 for the solver executable, found or fetched at a pinned revision.

The aggregate library target requires the Gurobi SDK even when the selected
runtime mode uses no Gurobi model.

The finder reads the version from `gurobi_c.h`, derives the core library name
(for example, `gurobi120` or `gurobi130`), and automatically searches for both
libraries in that SDK's `lib/` directory. It rejects versioned core libraries that
disagree with those headers. Explicit `GUROBI_LIBRARY` and `GUROBI_CXX_LIBRARY`
paths override the search; callers must keep those overrides consistent.
Set `GUROBI_HOME` to select an SDK; clear cached `GUROBI_*` paths when changing it.
Earlier discovery checks covered Gurobi 13 library names; it has not been compiled
or run here.

## Build and validation

Run from the repository root:

```sh
cmake -S . -B build
cmake --build build --target greedy-hps-moco
./build/greedy-hps-moco --help
```

`GHPS_BUILD_APPS` controls the solver executable only. CMake does not build,
register, or run tests. The MOKP experiment harness lives in `tests/mokp/`, moved
from `ignore/tests/`. Its reference generation, solver comparisons, and CSV output
are unchanged.

Compile the harness manually when needed. This example uses Gurobi 12 and the
default FetchContent dependency paths; set `GUROBI_HOME` to the SDK directory:

```sh
c++ -std=c++20 -O2 -fopenmp \
  -Iinclude -Ibuild/_deps/mooutils-src/include -Ibuild/_deps/cli11-src/include \
  -I"$GUROBI_HOME/include" tests/mokp/main.cpp \
  -L"$GUROBI_HOME/lib" -Wl,-rpath,"$GUROBI_HOME/lib" \
  -lgurobi_c++ -lgurobi120 -o build/greedy-hps-tests
./build/greedy-hps-tests --help
```

Adjust include paths for installed dependencies and the core library name for a
different SDK. The harness requires a supplied nondominated set, retains its
`--seed` option, and appends `results_*.csv` in its working directory. Gurobi solves
require a runtime license. `greedy-hps-tests` is a manually built executable, not
a CMake or install target.

## Install and consume the library

```sh
cmake -S . -B build-library -DGHPS_BUILD_APPS=OFF
cmake --install build-library --prefix /path/to/install
```

Consumers use `find_package(moco-hps CONFIG REQUIRED)` and link `moco::ghps`.
Set `CMAKE_PREFIX_PATH` to the installation prefix and `GUROBI_HOME` to the SDK.
The target supplies C++20, mooutils, OpenMP, and Gurobi. CLI dependencies and
project warning flags stay on executable targets.

Dataset downloading is optional: configure with `-DGHPS_FETCH_INSTANCES=ON`.
Instances then live under `build/_deps/mobkp_instances-src/instances/`. Read the
archived scripts before using them; `ignore/scripts/run_modes.sh` submits Slurm
jobs.

## Usage

```sh
./build/greedy-hps-moco --input-file=instance.in --mode=0 --J=10 --K=2 \
  --time-limit=60 --n-threads=2 --detailed-output=true
```

| Option | Meaning | Default |
| --- | --- | --- |
| `--input-file` | Instance path; required | — |
| `--mode` | Backend numbered 0–4 | `4` |
| `--J` | Maximum positive-gain selections | `5` |
| `--K` | Maximum reference bounds per scalarization; 1–30 | `3` |
| `--time-limit` | Positive time budget in seconds for optimization | `3600` |
| `--n-threads` | OpenMP subproblem workers; one Gurobi thread each | `1` |
| `--detailed-output` | Boolean labeled report | `false` |

Simulation requires a nonempty supplied front. The CLI rejects J greater than its
supplied size; the library's `GHPS_SIM` clamps oversized J. Simulation currently
has no time-budget parameter. Library optimization keeps the signature
`GHPS(instance, K, J, time_limit, n_threads, mode)`.

Timers use a monotonic clock and share deadlines across model construction,
enumeration, and optimization. Checks during construction are cooperative; a
single SDK call can return after the deadline.

## Input format

The parser reads whitespace-separated integers: `N M`, capacity, and N rows of
weight followed by M objective values. An optional final section contains a
nondominated-solution count followed by that many rows of M coordinates.
Comments are not part of the input syntax.

A complete two-item example is:

```text
2 2
1
1 4 2
1 2 4
2
4 2
2 4
```

The two rows after capacity describe items. The following `2` is the supplied
front size; `(4,2)` and `(2,4)` are its points. With J=1, simulation and BB return
one point, `HV(Nset)=12`, and `HV(Solutions)=8`; either point is a tied optimum.
Simulation writes a sequence file, so create its output directory first:

```sh
mkdir -p output
./build/greedy-hps-moco --input-file=instance.in --mode=4 --J=1 \
  --detailed-output=true
```

Coordinates and capacity are nonnegative int64 values; N must be positive and
M must be between 2 and 30. Objective sums must fit int64. These checks do not
guarantee that every accepted instance's hypervolume fits signed `__int128` or
that every integer coefficient is represented exactly in a floating-point model.

Optional benchmark data is available from
[the instance repository](https://github.com/gaplopes/mobkp-instances).

## Library structure and contracts

- `core/types.hpp` defines reference points, representation solutions, solver
  states, and 128-bit arithmetic helpers. `core/hypervolume.hpp` wraps mooutils
  and supports maximization and minimization reference regions.
- `core/nondominated_data.hpp` holds shared immutable front coordinates, direction,
  membership, and cached metrics, independently of knapsack data.
- `mokp/instance.hpp` owns one immutable item table in original input row order.
  `sorted_items[d]` contains row indices ordered by the ratio for objective d;
  `item_positions` maps rows to objective-0 search positions. BB traverses that
  sorted order and returns binary item flags in original row order.
- The direct constructor accepts N, M, capacity, Items, and front coordinates.
  The older constructor accepting separate weights and objective rows remains
  available and rejects disagreements with Items. Redundant public `weights`
  and `values` matrices are no longer stored; read `items[i].weight` and
  `items[i].values` instead.
- `search_region/local_bounds.hpp` maintains bounds; `subgraphs.hpp` constructs
  their adjacency and enumerates connected subsets with RwD. HPS reuses the
  adjacency across subset sizes and rebuilds it after bounds change.
- `hps/algorithm.hpp` coordinates optimization; `simulate.hpp` selects supplied
  candidates; `statistics.hpp` reports results. Instance and statistics snapshots
  share owned immutable front coordinates, membership, direction, and cached front
  metrics. Results survive the source instance. Metrics for selected solutions
  reflect the current public solution data.

Direction belongs to the instance and defaults to maximization. The MOKP
constructors and `MOKP::from_stream(stream, is_maximization)` accept it once;
Statistics and `GHPS_SIM(instance, J)` inherit it. Their former direction arguments
are removed, and `Statistics::is_maximization` is read-only. Both front volumes
are cached in that direction. The reference remains all-zero, so a minimization
front with nonnegative coordinates has zero origin-based hypervolume. Knapsack
optimization backends remain maximization-only and reject minimization instances.

Header names changed from `core/structs.hpp`, `core/hvi.hpp`, and
`search_region/llb.hpp` to the names above. Update consumer includes; no forwarding
headers are installed. `RefPoint`, `Solution`, and `LLB` type names are unchanged.

## Scores and output

`MOKPSolution` contains selected item flags, exact integer coordinates and weight,
termination status, backend status, node count, and a message. Its Gurobi
`hypervolume` is the **rounded model objective**: MILP reads ObjVal, MIQP rescales
ObjVal, and log-exp exponentiates ObjVal. BB computes its scalarization in integer
arithmetic. Floating-point model and extraction limits remain relevant.

`Solution.hypervolume` records the selected scalarization score.
`Solution.real_hypervolume` in optimization records the point's actual marginal
contribution to the complete representation; these are different quantities.
Simulation currently leaves that second metadata field zero. A separate exact
scalarization result field is not implemented in this change.

Statistics calculate representation hypervolume from selected integer coordinates.
Detailed output labels counts, timings, status, solutions, origin-based volumes,
and nadir-based volumes. Compact output still contains a status line and progress
output; it is not a CSV format. Per-subproblem tuples mean
`(reference_subset_size, node_count, elapsed_seconds)`.

`save_sequence(path)` writes the existing `.seq` format and checks open/write
errors. The CLI chooses `output/<instance>_<J>.seq`, preferring an existing
`output/<M>D/` directory. Reruns overwrite the same sequence file.

Gurobi's `optimal` is relative to its model tolerances. GHPS accepts completed
scalarizations and reports interrupted or failed computations with a nonzero CLI
exit status. Partial backend incumbents remain available for inspection. Tied
selections can differ between solvers or worker schedules.

The MOKP experiment harness compares rounded model scores with simulation;
its `hv_real` column currently repeats the reported model score. Large-integer
rounding differences can therefore produce disagreement rows. This remains a
separate unresolved review finding.

## Algorithm references

- [Klamroth, Lacour, and Vanderpooten: search-region representation](https://arxiv.org/abs/1502.06111).
- [Komusiewicz and Sommer: connected induced subgraph enumeration](https://www.uni-marburg.de/de/fb12/arbeitsgruppen/algorith/paper/enucon.pdf/@@download/file/enucon.pdf).

## License

MIT; see [LICENSE](LICENSE).
