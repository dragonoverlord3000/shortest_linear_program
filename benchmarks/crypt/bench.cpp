#include "bench.hpp"
#include "../io.hpp"

// the library
#include <slp/algorithm.hpp>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <chrono>
#include <exception>
#include <iostream>
#include <omp.h>
#include <string>

BenchResult run_crypt_benchmark(const Config &cfg) {
    auto t0_outer = std::chrono::steady_clock::now();

    std::vector<std::string> scheme_files =
        io::collect_recursive("benchmarks/crypt/dataset/");

    slp::Options options;
    options.alpha = cfg.potential_alpha;
    options.strategy = cfg.search_method;
    options.seed = cfg.seed;
    options.nearest = cfg.nearest;
    options.optimization_strategy = cfg.optimization_strategy;
    options.num_optimization_iters = cfg.num_optimization_iters;
    options.use_postprocess = cfg.use_postprocess;
    options.use_preprocess = cfg.use_preprocess;
    options.debug = cfg.debug;
    options.timelimit = cfg.timelimit;

    json root;
    root["config"] = {{"seed", cfg.seed},
                      {"num_basis_change", cfg.num_basis_change},
                      {"potential_alpha", cfg.potential_alpha},
                      {"potential_bit_p", cfg.potential_bit_p}};

    std::vector<json> scheme_results(scheme_files.size());
    std::exception_ptr first_exception = nullptr;

#pragma omp parallel for schedule(dynamic)
    for (std::ptrdiff_t i = 0;
         i < static_cast<std::ptrdiff_t>(scheme_files.size()); ++i) {
        try {
            const std::string &scheme_name = scheme_files[i];

            auto t0_inner = std::chrono::steady_clock::now();

            json scheme_json;
            scheme_json["scheme_name"] = scheme_name;
            scheme_json["results"] = json::object();

            auto [m, n, matrix] = io::parse_one_file_crypt(scheme_name);

            slp::Z2Matrix G(matrix, m, n);

            // Local copy in case slp::gf2::run mutates options internally.
            slp::Options local_options = options;

            slp::Result result = slp::gf2::run(G, local_options);

            auto t1_inner = std::chrono::steady_clock::now();

            if (cfg.verbose) {
#pragma omp critical(cout)
                {
                    std::cout << "scheme name: " << scheme_name << '\n';
                    std::cout
                        << "time: "
                        << std::chrono::duration_cast<std::chrono::nanoseconds>(
                               t1_inner - t0_inner)
                               .count()
                        << " ns\n";
                    std::cout << "best add for scheme [" << scheme_name
                              << "]: " << result.additions_after << '\n';
                }
            }

            json method_json = json::array();
            for (const auto &[a, b] : result.method.additions) {
                method_json.push_back({a, b});
            }

            auto duration_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(t1_inner -
                                                                      t0_inner)
                    .count();

            scheme_json["results"] = {{"best_add", result.additions_after},
                                      {"best_method", method_json},
                                      {"duration_ms", duration_ms}};

            scheme_results[i] = std::move(scheme_json);
        } catch (...) {
#pragma omp critical(exception_capture)
            {
                if (!first_exception) {
                    first_exception = std::current_exception();
                }
            }
        }
    }

    if (first_exception) {
        std::rethrow_exception(first_exception);
    }

    root["schemes"] = json::array();
    for (auto &scheme_json : scheme_results) {
        root["schemes"].push_back(std::move(scheme_json));
    }

    auto t1_outer = std::chrono::steady_clock::now();
    auto wall_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                            t1_outer - t0_outer)
                            .count();

    BenchResult bench_result;
    bench_result.name = "crypt";
    bench_result.wall_time_ms = wall_time_ms;
    bench_result.instances = scheme_files.size();
    bench_result.details = std::move(root);
    return bench_result;
}
