#include "iostream"
#include <thread>

#include "3x3_matmul/bench.hpp"
#include "types.hpp"

// compiler + c++ flags the benchmark is run with
#ifndef SLP_COMPILER
#define SLP_COMPILER "unknown"
#endif

#ifndef SLP_CXXFLAGS
#define SLP_CXXFLAGS "unknown"
#endif

int main() {
    Config cfg;

    std::vector<BenchResult> results;

    std::cout << "running 3x3 matmul benchmark..." << std::endl;
    results.push_back(run_3x3_matmul_benchmark(cfg));


    json j;
    j["config"] = {
        {"search_method", cfg.search_method},
        {"output", cfg.output},
        {"threads", cfg.threads},
        {"seed", cfg.seed},
        {"3x3_matmul_details", {cfg.num_basis_change}}
    };
    j["build"] = {
        {"compiler", SLP_COMPILER},
        {"cxxflags", SLP_CXXFLAGS}
    };
    j["machine"] = {
        {"logical_cpus", std::thread::hardware_concurrency()}
    };

    // TODO
    j["results"];

    return 0;
}
