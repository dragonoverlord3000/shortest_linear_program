#include <fstream>
#include <iostream>
#include <thread>

#include <argparse/argparse.hpp>

#include "3x3_matmul/bench.hpp"
#include "crypt/bench.hpp"
#include "struct_matmul/bench.hpp"

#include "slp/types.hpp"
#include "types.hpp"

// compiler + c++ flags the benchmark is run with
#ifndef SLP_COMPILER
#define SLP_COMPILER "unknown"
#endif

#ifndef SLP_CXXFLAGS
#define SLP_CXXFLAGS "unknown"
#endif

void add_arguments(argparse::ArgumentParser &program, Config &cfg) {
    // general
    program.add_argument("--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true)
        .store_into(cfg.verbose);
    program.add_argument("--debug")
        .help("increase output verbosity by a lot")
        .default_value(false)
        .implicit_value(true)
        .store_into(cfg.debug);

    program.add_argument("--benchmarks")
        .help("which benchmarks to run")
        .choices(std::string{"3x3_matmul"}, std::string{"crypt"},
                 std::string{"struct_matmul"})
        .nargs(argparse::nargs_pattern::at_least_one)
        .store_into(cfg.benchmarks);
    program.add_argument("--output")
        .help("set the output json path")
        .default_value(std::string("build/benchmarks/full.json"))
        .nargs(1)
        .store_into(cfg.output);
    program.add_argument("--seed")
        .help("specify the seed to use, same seed is shared everywhere")
        .default_value(uint64_t{628318})
        .nargs(1)
        .scan<'u', uint64_t>()
        .store_into(cfg.seed);
    program.add_argument("--max_seconds")
        .help("maximum time spent on any of the benchmarks")
        .default_value(std::numeric_limits<uint64_t>::max())
        .nargs(1)
        .scan<'u', uint64_t>()
        .store_into(cfg.max_seconds);

    program.add_argument("--search_method")
        .help("set the heuristic to use, can be one of")
        .default_value(std::string("greedy_potential"))
        .choices("greedy_potential", "backtrack_potential", "BP", "RNBP", "A1",
                 "A2", "paar1")
        .nargs(1);
    program.add_argument("--timelimit")
        .help("the timelimit for each matrix")
        .default_value(0.10) // 0.10 seconds
        .nargs(1)
        .scan<'g', double>()
        .store_into(cfg.timelimit);

    program.add_argument("--optimization_strategy")
        .help("use 'framework', 'single_shot', or 'repeat_random' to minimize "
              "additions")
        .default_value(std::string("framework"))
        .choices("framework", "single_shot", "repeat_random")
        .nargs(1);
    program.add_argument("--num_optimization_iters")
        .help("how many iterations of framework to run")
        .default_value(size_t{std::numeric_limits<size_t>::max()})
        .nargs(1)
        .scan<'u', size_t>()
        .store_into(cfg.num_optimization_iters);

    // specific for potential methods
    program.add_argument("--potential_alpha")
        .help("how much weight to put on potential")
        .default_value(0.2)
        .nargs(1)
        .scan<'g', double>()
        .store_into(cfg.potential_alpha);

    // specific for A1, A2 methods
    program.add_argument("--ax_nearest")
        .help("how relaxed to have filtering step in A1, A2 heuristics")
        .default_value(size_t{0})
        .nargs(1)
        .scan<'u', size_t>()
        .store_into(cfg.nearest);

    // specific for 3x3 matmul benchmark
    program.add_argument("--num_basis_change")
        .help("how many basis change matrices to use per 3x3 matmul scheme")
        .default_value(uint64_t{1})
        .nargs(1)
        .scan<'u', uint64_t>()
        .store_into(cfg.num_basis_change);
    program.add_argument("--potential_bit_p")
        .help("basis change matrices are sampled from GL bernouli random "
              "matrices with bit_p parameters")
        .default_value(0.15)
        .nargs(1)
        .scan<'g', double>()
        .store_into(cfg.potential_bit_p);
    program.add_argument("--specific_type")
        .help("whether to use all types (empty), or only 'W', 'U' or 'V'")
        .default_value(std::string(""))
        .nargs(1)
        .store_into(cfg.specific_type);

    // specific for crypt
    //
    // ...

    // hardware
    program.add_argument("--threads")
        .help("number of threads the job is run with")
        .default_value(uint32_t{1})
        .nargs(1)
        .scan<'u', uint32_t>()
        .store_into(cfg.threads);
}

// for those arguments that cannot be filled in directly
void fill_cfg(argparse::ArgumentParser &program, Config &cfg) {
    std::string search_method = program.get<std::string>("--search_method");
    std::string optimization_strategy =
        program.get<std::string>("--optimization_strategy");

    if (search_method == "greedy_potential") {
        cfg.search_method = slp::SearchStrategy::GreedyPotential;
    } else if (search_method == "backtrack_potential") {
        cfg.search_method = slp::SearchStrategy::BacktrackingPotential;
    } else if (search_method == "BP") {
        cfg.search_method = slp::SearchStrategy::BP;
    } else if (search_method == "paar1") {
        cfg.search_method = slp::SearchStrategy::Paar1;
    } else if (search_method == "RNBP") {
        cfg.search_method = slp::SearchStrategy::RNBP;
    } else if (search_method == "A1") {
        cfg.search_method = slp::SearchStrategy::A1;
    } else if (search_method == "A2") {
        cfg.search_method = slp::SearchStrategy::A2;
    } else {
        throw std::invalid_argument(
            "received invalid argument for search method");
    }

    if (optimization_strategy == "framework") {
        cfg.optimization_strategy = slp::OptimizationStrategy::Framework;
    } else if (optimization_strategy == "single_shot") {
        cfg.optimization_strategy = slp::OptimizationStrategy::SingleShot;
    } else if (optimization_strategy == "repeat_random") {
        cfg.optimization_strategy = slp::OptimizationStrategy::RepeatRandom;
    } else {
        throw std::invalid_argument(
            "received invalid argument for optimization strategy");
    }
}

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("benchmarker");
    Config cfg;
    add_arguments(program, cfg);
    try {
        program.parse_args(argc, argv);
        fill_cfg(program, cfg);
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    // some main stats
    std::cout << "running the benchmarks..." << std::endl;
    std::vector<BenchResult> results;
    for (const std::string &benchmark : cfg.benchmarks) {
        if (benchmark == "3x3_matmul") {
            std::cout << "running 3x3 matmul benchmark..." << std::endl;
            results.push_back(run_3x3_matmul_benchmark(cfg));
        } else if (benchmark == "crypt") {
            std::cout << "running crypt benchmarks..." << std::endl;
            results.push_back(run_crypt_benchmark(cfg));
        } else if (benchmark == "struct_matmul") {
            std::cout
                << "running structured matrix multiplication benchmarks..."
                << std::endl;
            results.push_back(run_struct_benchmark(cfg));
        } else {
            throw std::runtime_error("unknown benchmark: " + benchmark);
        }
    }

    // record all the benchmark outputs
    json j;
    j["config"] = {
        {"search_method", program.get<std::string>("--search_method")},
        {"optimization_strategy",
         program.get<std::string>("--optimization_strategy")},
        {"method_details",
         {{"potential_alpha", cfg.potential_alpha},
          {"nearest", cfg.nearest},
          {"num_optimization_iters", cfg.num_optimization_iters},
          {"timelimit", cfg.timelimit}}},
        {"output", cfg.output},
        {"threads", cfg.threads},
        {"seed", cfg.seed},
        {"benchmarks", cfg.benchmarks},
        {"3x3_matmul_details",
         {{"num_basis_change", cfg.num_basis_change},
          {"potential_bit_p", cfg.potential_bit_p},
          {"specific_type", cfg.specific_type}}}};
    j["build"] = {{"compiler", SLP_COMPILER}, {"cxxflags", SLP_CXXFLAGS}};
    j["machine"] = {{"logical_cpus", std::thread::hardware_concurrency()}};
    j["results"] = json::array();
    for (BenchResult &result : results) {
        json j_result;
        j_result["name"] = result.name;
        j_result["instances"] = result.instances;
        j_result["wall_time_ms"] = result.wall_time_ms;
        j_result["details"] = result.details;
        j["results"].push_back(j_result);
    }

    std::ofstream out(cfg.output);
    if (!out) {
        std::cerr << "failed to open output file: " << cfg.output << std::endl;
        return 1;
    }
    out << j.dump(2) << std::endl;

    return 0;
}
