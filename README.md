This directory is for all the production-ready code. 
Much of it is copied from the relevant parts of the `exploration` folder.

# Algorithms
## Greedy Potential
The repository includes a fast implementation of the greedy potential method for reducing the number of additions in the matrix-vector product $Gx$ for either 
$G \in \{-1, 0, 1\}^{m \times n}$ or $G \in \mathbb{Z}_2$ described in my MSc Thesis.

Inspiration: (https://eprint.iacr.org/2024/2063.pdf)


## Boyar-Peralta 
Inspiration: (https://eprint.iacr.org/2025/1493.pdf)


# Quick Start
First run 
```bash
make
```
to build the algorithm. To then use it in your own work, do:
TODO ???


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








