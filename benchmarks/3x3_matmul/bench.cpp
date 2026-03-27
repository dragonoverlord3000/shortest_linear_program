#include "bench.hpp"
#include "io.hpp"
#include "matrix.hpp"

// the library
#include <slp/algorithm.hpp>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <unordered_set>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <chrono>
#include <limits>

// small utility function for when bit_p is very low, so that we get unique
// basis change matrices
std::size_t array_hasher(std::array<uint16_t, 9> const &vec) {
    std::size_t seed = vec.size();
    for (auto &i : vec) {
        seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

BenchResult run_3x3_matmul_benchmark(const Config &cfg) {
    auto t0_outer = std::chrono::steady_clock::now();
    std::vector<std::string> scheme_files =
        io::collect_recursive("benchmarks/3x3_matmul/schemes-tab/schemes/"); // assumes running from root
    std::vector<std::string> types = {"U", "V", "W"};
    if (!cfg.specific_type.empty()) types = {cfg.specific_type};

    slp::Options options;
    options.alpha = cfg.potential_alpha;
    options.strategy = cfg.search_method;

    std::mt19937 rng(cfg.seed);
    std::vector<std::array<uint16_t, 9>> Bs(cfg.num_basis_change);

    // let Bs[0] be the identity matrix
    for (int i = 0; i < 9; i++)
        Bs[0][i] = 1 << i;
    std::unordered_set<size_t> repeat_counter;
    repeat_counter.insert(array_hasher(Bs[0]));
    for (std::size_t i = 1; i < cfg.num_basis_change; i++) {
        do {
            matrix::fill_random_GL9(rng, Bs[i], cfg.potential_bit_p);
        } while (repeat_counter.count(array_hasher(Bs[i])));
        repeat_counter.insert(array_hasher(Bs[i]));
    }

    json root;
    root["config"] = {{"seed", cfg.seed},
                      {"num_basis_change", cfg.num_basis_change},
                      {"potential_alpha", cfg.potential_alpha},
                      {"potential_bit_p", cfg.potential_bit_p}};
    root["schemes"] = json::array();

    uint64_t instances = 0;
    for (std::size_t scheme_idx = 0; scheme_idx < scheme_files.size();
         scheme_idx++) {
        json scheme_json;
        scheme_json["scheme_name"] = scheme_files[scheme_idx];
        scheme_json["results"] = json::object();

        for (const std::string& type : types) {
            auto t0_inner = std::chrono::steady_clock::now();
            std::size_t m = type == "W" ? 9 : 23;
            std::size_t n = type == "W" ? 23 : 9;
            std::vector<uint64_t> G =
                io::parse_one_file(scheme_files[scheme_idx], type);

            std::size_t best_add = std::numeric_limits<std::size_t>::max();
            slp::AdditionMethod best_method;
            std::size_t best_basis_change_idx = 0;

            for (std::size_t rand_idx = 0; rand_idx < cfg.num_basis_change;
                 rand_idx++) {
                instances++;
                std::vector<uint64_t> BG;
                BG.reserve(G.size());
                if (type == "W")
                    matrix::_change_basis_W(G, Bs[rand_idx], BG);
                else
                    matrix::_change_basis_UV(G, Bs[rand_idx], BG);

                slp::Z2Matrix _G(BG, m, n);
                slp::Result result = slp::gf2::run(_G, options);

                if (result.additions_after < best_add) {
                    best_add = result.additions_after;
                    best_method = result.method;
                    best_basis_change_idx = rand_idx;
                }
            }
            auto t1_inner = std::chrono::steady_clock::now();

            if (cfg.verbose) {
                std::cout << "scheme idx: " << scheme_idx << std::endl;
                std::cout << "type: " << type << std::endl;
                std::cout << "time: " << static_cast<std::chrono::nanoseconds>(t1_inner - t0_inner) << std::endl;
                std::cout << "best add for scheme [" << scheme_idx << "]: " << best_add << std::endl;
            }

            json method_json = json::array();
            for (const auto &[a, b] : best_method.additions) {
                method_json.push_back({a, b});
            }
            
            auto duration_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(t1_inner - t0_inner).count();
            scheme_json["results"][type] = {
                {"best_add", best_add},
                {"best_basis_change_idx", best_basis_change_idx},
                {"best_basis_change", Bs[best_basis_change_idx]},
                {"best_method", method_json},
                {"duration_ms", duration_ms}};
        }

        root["schemes"].push_back(std::move(scheme_json));
    }

    auto t1_outer = std::chrono::steady_clock::now();
    auto wall_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1_outer - t0_outer).count();

    BenchResult bench_result;
    bench_result.name = "3x3_matmul";
    bench_result.wall_time_ms = wall_time_ms; 
    bench_result.instances = instances;
    bench_result.details = std::move(root);
    return bench_result;
}
