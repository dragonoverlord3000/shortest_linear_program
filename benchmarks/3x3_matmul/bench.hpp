#pragma once
#include "../types.hpp"

std::size_t array_hasher(std::array<uint16_t, 9> const &vec);
BenchResult run_3x3_matmul_benchmark(const Config& cfg);
