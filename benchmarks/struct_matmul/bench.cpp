#include "bench.hpp"
#include "../io.hpp"

// the library
#include <slp/algorithm.hpp>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <chrono>
#include <iostream>
#include <random>
#include <string>

BenchResult run_struct_benchmark(const Config &cfg) {
    auto t0_outer = std::chrono::steady_clock::now();
    std::vector<std::string> scheme_files = io::collect_recursive(
        "benchmarks/struct_matmul/dataset/"); // assumes running from root

    slp::Options options;
    options.alpha = cfg.potential_alpha;
    options.strategy = cfg.search_method;
    options.seed = cfg.seed;
    options.nearest = cfg.nearest;
    options.optimization_strategy = cfg.optimization_strategy;
    options.num_optimization_iters = cfg.num_optimization_iters;

    options.debug = cfg.debug;
    options.timelimit = cfg.timelimit;

    std::mt19937 rng(cfg.seed);

    json root;
    root["config"] = {{"seed", cfg.seed},
                      {"num_basis_change", cfg.num_basis_change},
                      {"potential_alpha", cfg.potential_alpha},
                      {"potential_bit_p", cfg.potential_bit_p}};
    root["schemes"] = json::array();

    uint64_t instances = 0;
    for (const std::string &scheme_name : scheme_files) {
        if (scheme_name.size() < 5 ||
            scheme_name.substr(scheme_name.size() - 4) != ".txt" ||
            (scheme_name[scheme_name.size() - 5] != 'A' &&
             scheme_name[scheme_name.size() - 5] != 'B' &&
             scheme_name[scheme_name.size() - 5] != 'M'))
            continue;
        auto t0_inner = std::chrono::steady_clock::now();
        instances++;

        json scheme_json;
        scheme_json["scheme_name"] = scheme_name;
        scheme_json["results"] = json::object();
        auto [m, n, matrix] = io::parse_one_file_struct(scheme_name);
        slp::Z2Matrix G(matrix, m, n);
        slp::Result result = slp::gf2::run(G, options);

        auto t1_inner = std::chrono::steady_clock::now();
        if (cfg.verbose) {
            std::cout << "scheme name: " << scheme_name << std::endl;
            std::cout << "time: "
                      << static_cast<std::chrono::nanoseconds>(t1_inner -
                                                               t0_inner)
                      << std::endl;
            std::cout << "best add for scheme [" << scheme_name
                      << "]: " << result.additions_after << std::endl;
        }

        json method_json = json::array();
        for (const auto &[a, b] : result.method.additions)
            method_json.push_back({a, b});

        auto duration_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(t1_inner -
                                                                  t0_inner)
                .count();
        scheme_json["results"] = {{"best_add", result.additions_after},
                                  {"best_method", method_json},
                                  {"duration_ms", duration_ms}};

        root["schemes"].push_back(std::move(scheme_json));
    }

    auto t1_outer = std::chrono::steady_clock::now();
    auto wall_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                            t1_outer - t0_outer)
                            .count();

    BenchResult bench_result;
    bench_result.name = "3x3_matmul";
    bench_result.wall_time_ms = wall_time_ms;
    bench_result.instances = instances;
    bench_result.details = std::move(root);
    return bench_result;
}
