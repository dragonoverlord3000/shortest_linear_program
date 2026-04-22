#pragma once

#include "slp/types.hpp"

#include <stdint.h>
#include <string>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

struct Config {
    std::vector<std::string> benchmarks = {"3x3_matmul"};
    slp::SearchStrategy search_method = slp::SearchStrategy::GreedyPotential;
    slp::OptimizationStrategy optimization_strategy = slp::OptimizationStrategy::Framework;
    size_t num_optimization_iters = 32;
    double prob_framework_include = 0.4;

    // specifically for potential method
    double potential_alpha = 0.2;

    // specifically for 3x3 matmul benchmark
    uint64_t num_basis_change = 1;
    double potential_bit_p = 0.25;
    std::string specific_type =
        ""; // empty means all. Other possible choices are "W", "U", "V"

    // specifically for crypt
    // ...

    // for A1, A2 algorithms
    size_t nearest = 0;

    std::string output = "build/benchmarks/full.json";
    uint64_t seed = 628318;
    uint64_t max_seconds = std::numeric_limits<uint64_t>::max();
    bool verbose = true;
    uint32_t threads = 1;
};

struct BenchResult {
    std::string name;
    uint64_t wall_time_ms = 0;
    uint64_t instances = 0;
    json details;
};
