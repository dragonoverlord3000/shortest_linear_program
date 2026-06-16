# Shortest Linear Program 
`C++` library of methods for optimizing instances of the Shortest Linear Program (SLP) problem. Specifically, given a matrix 

$$ 
G \in \mathbb{F}_2^{m \times n}.
$$

We find low-addition methods for computing the matrix-vector product $Gx$ over $\mathbb{F}_2$ that work for all 
vectors $x \in \mathbb{F}_2^{n}$.

The code is part of my MSc thesis, see `thesis.pdf`.

## Algorithms

The following search methods are currently available for binary matrices over \(\mathbb{F}_2\):

| CLI name | Description |
|---|---|
| `greedy_potential` | Greedy potential-based heuristic. |
| `backtrack_potential` | Backtracking variant of the potential method. |
| `BP` | Boyar-Peralta inspired heuristic. |
| `RNBP` | Randomized variant of the Boyar-Peralta inspired heuristic. |
| `A1` | Ax-style heuristic, variant A1. |
| `A2` | Ax-style heuristic, variant A2. |
| `paar1` | Paar-style common-subexpression heuristic. |

For ternary matrices over $\{-1,0,1\}$, only the greedy potential method is currently supported.


## Build
The default build does not include ORTools, making it so that the mixed integer programming is unavailable.
To run the default build, simply do:

```bash
make
```


### Notes
On LSF systems, load a recent GCC before building:

```bash
module load gcc/15.2.0-binutils-2.45
make
```

If you want to build with ORTools, and it is installed somewhere not on the default path, then 
use:
```bash
make USE_MATHOPT=1 ORTOOLS_DIR=/path/to/or-tools
```

## CLI usage
To run the CLI code, do:
```bash
./build/bin/slp [options] < input.txt
```
where `input.txt` is on the form:
```text
m n
g_{0,0} g_{0,1} ... g_{0,n-1}
g_{1,0} g_{1,1} ... g_{1,n-1}
...
g_{n-1,0} g_{n-1,1} ... g_{n-1,n-1}
```
with all entries outside $m$ and $n$ being either $0$ or $1$. For example `input.txt` might look like:
```text
2 4
1 1 1 0
0 1 1 1
```

### Ternary
To run the single ternary implementation run:
```bash
./build/bin/slp --ternary --optimization_strategy single_shot --search_method greedy_potential < input.txt
```
With the entries of the matrix in `input.txt` being $0$, $1$ or $-1$.

Currently ternary only supports 
```bash
--optimization_strategy single_shot
--search_method greedy_potential
```

### CLI Options
| Option                     |                                                                              Values |                       Default | Description                                                                                                                                                                      |
| -------------------------- | ----------------------------------------------------------------------------------: | ----------------------------: | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--verbose`                |                                                                                flag |                       `false` | Increase output verbosity.                                                                                                                                                       |
| `--debug`                  |                                                                                flag |                       `false` | Print much more debugging information. Recommended only when debugging.                                                                                                          |
| `--skip_first`             |                                                                                flag |                       `false` | Skip the first integer in the input stream before reading the matrix.                                                                                                            |
| `--seed`                   |                                                                    unsigned integer |                      `628318` | Random seed shared across randomized parts of the solver.                                                                                                                        |
| `--no-preprocess`          |                                                                                flag |         preprocessing enabled | Disable preprocessing.                                                                                                                                                           |
| `--no-postprocess`         |                                                                                flag |        postprocessing enabled | Disable postprocessing.                                                                                                                                                          |
| `--search_method`          | `greedy_potential`, `backtrack_potential`, `BP`, `RNBP`, `A1`, `A2`, `paar1`, `MIP` |            `greedy_potential` | Search heuristic to use.                                                                                                                                                         |
| `--reachable_strategy`     |                                `backtracking_sparsity_aware`, `brute_force`, `mitm` | `backtracking_sparsity_aware` | Reachability strategy used by BP-inspired heuristics.                                                                                                                            |
| `--timelimit`              |                                                              floating-point seconds |                        `60.0` | Time limit for solving the matrix.                                                                                                                                               |
| `--optimization_strategy`  |                                         `framework`, `single_shot`, `repeat_random` |                   `framework` | Outer optimization strategy used to minimize additions.                                                                                                                          |
| `--num_optimization_iters` |                                                                    unsigned integer |                  max `size_t` | Number of framework iterations to run.                                                                                                                                           |
| `--potential_alpha`        |                                                               floating-point number |                         `0.2` | Weight used by the potential heuristic.                                                                                                                                          |
| `--ax_nearest`             |                                                                    unsigned integer |                           `0` | Relaxation parameter used by the A1/A2 filtering step.                                                                                                                           |
| `--ternary`                |                                                                                flag |                       `false` | Use ternary matrix mode over ${-1,0,1}$ instead of binary $\mathbb{F}_2$ mode. Only supported with `--optimization_strategy single_shot` and `--search_method greedy_potential`. |

## Output Format
The initial ordered basis $B$ is
$$
B_0 = e_0, \quad B_1 = e_1, \quad \ldots, \quad B_{n-1} = e_{n-1}.
$$
Each ($0$-indexed) printed pair $i,j$ appends to $B$ the value of $B_i \oplus B_j$, namely
$$
B_{|B| + 1} \gets B_i \oplus B_j.
$$
The output indices specify which basis vector computes each output row of $Gx$.

We use `std::numeric_limits<size_t>::max()` to mark zero-rows.

## Library Usage
Include the public headers from `include/` and link against `build/lib/libslp.a`.

Example compile command:
```bash
g++ -std=c++23 -Iinclude -Ithird_party your_program.cpp \
    -Lbuild/lib -lslp \
    -o your_program
```

Use the examples provided in `examples/`, and the corresponding part of the `Makefile` as a guide.


## Benchmarks
To download a certain benchmark dataset, run
```
make download-x
```
with `x` being one of `3x3`, `crypt`, `struct`, or `bernoulli`. Note that the `3x3` dataset has some sort of 
rate limiting on their website, so either change your IP address before re-downloading, or just wait a long time.

To **build and run** the benchmarks do
```bash
make bench-full BENCH_ARGS="..."
```
Where possible benchmark arguments include:

| Argument                   | Type                                                                                |                                Default | Description                                                                  |
| -------------------------- | ----------------------------------------------------------------------------------- | -------------------------------------: | ---------------------------------------------------------------------------- |
| `--verbose`                | flag                                                                                |                                `false` | Increase output verbosity.                                                   |
| `--debug`                  | flag                                                                                |                                `false` | Increase output verbosity by a lot.                                          |
| `--benchmarks`             | one or more of `3x3_matmul`, `crypt`, `struct_matmul`, `bernoulli`                  |                                      — | Select which benchmarks to run.                                              |
| `--output`                 | path                                                                                |           `build/benchmarks/full.json` | Path to the output JSON file.                                                |
| `--seed`                   | unsigned integer                                                                    |                               `628318` | Random seed shared across all benchmarks.                                    |
| `--max_seconds`            | unsigned integer                                                                    | `std::numeric_limits<uint64_t>::max()` | Maximum time spent on any benchmark.                                         |
| `--search_method`          | `greedy_potential`, `backtrack_potential`, `BP`, `RNBP`, `A1`, `A2`, `paar1`, `MIP` |                     `greedy_potential` | Search heuristic to use.                                                     |
| `--reachable_strategy`     | `backtracking_sparsity_aware`, `brute_force`, `mitm`                                |          `backtracking_sparsity_aware` | Strategy for finding reachability in BP-inspired heuristics.                 |
| `--timelimit`              | floating-point number                                                               |                                 `0.10` | Time limit, in seconds, for each matrix.                                     |
| `--optimization_strategy`  | `framework`, `single_shot`, `repeat_random`                                         |                            `framework` | Optimization strategy used to minimize additions.                            |
| `--num_optimization_iters` | unsigned integer                                                                    |   `std::numeric_limits<size_t>::max()` | Number of framework optimization iterations to run.                          |
| `--no-preprocess`          | flag                                                                                |                                `false` | Disable preprocessing.                                                       |
| `--no-postprocess`         | flag                                                                                |                                `false` | Disable postprocessing.                                                      |
| `--potential_alpha`        | floating-point number                                                               |                                  `0.2` | Weight assigned to potential in the heuristic.                               |
| `--ax_nearest`             | unsigned integer                                                                    |                                    `0` | Relaxation parameter for the filtering step in the `A1` and `A2` heuristics. |
| `--num_basis_change`       | unsigned integer                                                                    |                                    `1` | Number of basis-change matrices to use per `3x3_matmul` scheme.              |
| `--potential_bit_p`        | floating-point number                                                               |                                 `0.15` | Bernoulli parameter used when sampling basis-change matrices.                |
| `--specific_type`          | `""`, `W`, `U`, `V`                                                                 |                                   `""` | Restrict to a specific type, or leave empty to use all types.                |
| `--threads`                | unsigned integer                                                                    |                                    `1` | Number of threads to run with.                                               |



# Repository Structure
```bash
.
├── apps/
├── benchmarks/
│   ├── 3x3_matmul/
│   ├── bernoulli/
│   ├── crypt/
│   └── struct_matmul/
├── examples/
├── extender/
├── include/
│   └── slp/
│       ├── boyar_peralta/
│       ├── framework/
│       ├── mip/
│       ├── paar/
│       ├── postprocess/
│       ├── potential/
│       ├── preprocess/
│       └── utils/
├── src/
│   └── slp/
│       ├── boyar_peralta/
│       ├── framework/
│       ├── mip/
│       ├── paar/
│       ├── postprocess/
│       ├── potential/
│       ├── preprocess/
│       └── utils/
└── third_party/
    ├── argparse/
    └── nlohmann/
```








