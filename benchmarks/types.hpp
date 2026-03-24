#pragma once

#include "slp/types.hpp"
#include <string>
#include <stdint.h>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

struct Config {
    slp::SearchStrategy search_method = slp::SearchStrategy::GreedyPotential;

    // specifically for potential method
    double potential_alpha = 0.2;
    double potential_bit_p = 0.25;

    // specifically for 3x3 matmul benchmark
    uint64_t num_basis_change = 1;

    // specifically for crypt
    // ...

    std::string output = "build/benchmarks/full.json";
    uint64_t seed = 628318;
    uint32_t threads = 1;
};

struct BenchResult {
    std::string name;
    uint64_t wall_time_ms = 0;
    uint64_t instances = 0;
    json details;
};
