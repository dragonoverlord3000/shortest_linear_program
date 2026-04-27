This directory is for all the production-ready code. 
Much of it is copied from the relevant parts of the `exploration` folder.


On LSF remember:
```
module load gcc/15.2.0-binutils-2.45
```


# Algorithms
## Greedy Potential
The repository includes a fast implementation of the greedy potential method for reducing the number of additions in the matrix-vector product $Gx$ for either 
$G \in \{-1, 0, 1\}^{m \times n}$ or $G \in \mathbb{Z}_2$ described in my MSc Thesis.

Inspiration: (https://eprint.iacr.org/2024/2063.pdf)


## Boyar-Peralta 
Inspiration: (https://eprint.iacr.org/2025/1493.pdf)


# Quick Start
First clone this repository or add it as a submodule to your own project:
```bash
git clone ...
```
To build the library, run 
```bash
make
```
Then to use it in your own project do:
TODO ???

## Benchmarks
To **build and run** the benchmarks do
```bash
make bench-full BENCH_ARGS="..."
```
Where possible benchmark arguments include:

| Argument | Type / Values | Default | Description |
|---|---|---:|---|
| `--verbose` | flag | `false` | Increase output verbosity. |
| `--benchmarks` | one or more of `3x3_matmul`, `crypt` | — | Select which benchmarks to run. |
| `--output` | path | `build/benchmarks/full.json` | Path to the output JSON file. |
| `--seed` | unsigned integer | `628318` | Random seed shared across all benchmarks. |
| `--search_method` | `greedy_potential`, `backtrack_potential` | `greedy_potential` | Search heuristic to use. |
| `--potential_alpha` | floating-point number | `0.2` | Weight assigned to potential in the heuristic. |
| `--num_basis_change` | unsigned integer | `1` | Number of basis-change matrices to use per `3x3_matmul` scheme. |
| `--potential_bit_p` | floating-point number | `0.25` | Bernoulli parameter used when sampling basis-change matrices. |
| `--specific_type` | `""`, `W`, `U`, `V` | `""` | Restrict to a specific type, or leave empty to use all types. |
| `--threads` | unsigned integer | `1` | Number of threads to run with. |



# Repository Structure
```bash
.
├── Makefile
├── README.md
├── compile_commands.json
├── include/
│   └── slp/
│       ├── algorithm.hpp
│       ├── types.hpp
│       └── potential/
│           └── internal.hpp
├── src/
│   └── slp/
│       ├── algorithm_gf2.cpp
│       ├── boyar_peralta/
│       └── potential/
│           └── gf2/
│               ├── backtrack.cpp
│               ├── core.cpp
│               └── greedy.cpp
├── examples/
│   └── use_gf2.cpp
├── benchmarks/
└── build/
    ├── examples/
    ├── lib/
    └── obj/
```








