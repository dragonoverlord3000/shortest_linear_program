#include <iostream>

#include <argparse/argparse.hpp>

#include "slp/algorithm.hpp"
#include "slp/potential/internal.hpp"
#include "slp/types.hpp"

namespace {
void read_binary_matrix(size_t &m, size_t &n, std::vector<uint64_t> &matrix) {
    std::cin >> m >> n;
    matrix.assign(n, 0ULL);

    for (size_t i = 0; i < m; i++) {
        uint64_t bit = 1ULL << i;
        for (size_t j = 0; j < n; j++) {
            uint64_t x;
            std::cin >> x;
            assert(x == 0 || x == 1);
            matrix[j] |= x * bit;
        }
    }
}

int get_ternary_naive_adds(std::vector<std::vector<int>> &matrix) {
    int ans = 0;
    for (std::vector<int> &row : matrix) {
        int cnt = 0;
        for (int v : row)
            cnt += v != 0;
        ans += std::max(0, cnt - 1);
    }
    return ans;
}

void read_ternary_matrix(size_t &m, size_t &n,
                         std::vector<std::vector<int>> &matrix) {
    std::cin >> m >> n;
    matrix.assign(m, std::vector<int>(n, 0));
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            std::cin >> matrix[i][j];
}

void validate_options(const slp::Options &options) {
    if (options.ternary) {
        if (options.optimization_strategy !=
                slp::OptimizationStrategy::SingleShot ||
            options.strategy != slp::SearchStrategy::GreedyPotential) {
            throw std::invalid_argument(
                "--ternary is only supported with "
                "--optimization-strategy single_shot and "
                "--search-method greedy_potential");
        }
    }
}

} // namespace

void add_arguments(argparse::ArgumentParser &program, slp::Options &options) {
    // general
    program.add_argument("--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true)
        .store_into(options.verbose);
    program.add_argument("--debug")
        .help("increase output verbosity by a lot, recommend having it off "
              "unless debugging")
        .default_value(false)
        .implicit_value(true)
        .store_into(options.debug);

    program.add_argument("--skip_first")
        .help("skip the first number when reading the matrix")
        .default_value(false)
        .implicit_value(true)
        .store_into(options.skip_first);

    program.add_argument("--seed")
        .help("specify the seed to use, same seed is shared everywhere")
        .default_value(uint64_t{628318})
        .nargs(1)
        .scan<'u', uint64_t>()
        .store_into(options.seed);

    program.add_argument("--no-preprocess")
        .help("disable preprocessing")
        .flag()
        .action([&](const auto &) { options.use_preprocess = false; });

    program.add_argument("--no-postprocess")
        .help("disable postprocessing")
        .flag()
        .action([&](const auto &) { options.use_postprocess = false; });

    program.add_argument("--search_method")
        .help("set the heuristic to use, can be one of")
        .default_value(std::string("greedy_potential"))
        .choices("greedy_potential", "backtrack_potential", "BP", "RNBP", "A1",
                 "A2", "paar1", "MIP")
        .nargs(1);
    program.add_argument("--reachable_strategy")
        .help("set the strategy for finding reachability in BP inspired "
              "heuristics")
        .default_value("backtracking_sparsity_aware")
        .choices("backtracking_sparsity_aware", "brute_force", "mitm")
        .nargs(1)
        .action([&](const auto &reachable_strategy) {
            if (reachable_strategy == "backtracking_sparsity_aware") {
                options.reachable_strategy =
                    slp::ReachableStrategy::BacktracingSparseAware;
            } else if (reachable_strategy == "brute_force") {
                options.reachable_strategy = slp::ReachableStrategy::BruteForce;
            } else if (reachable_strategy == "mitm") {
                options.reachable_strategy = slp::ReachableStrategy::MITM;
            }
        });

    program.add_argument("--timelimit")
        .help("the timelimit for solving the matrix")
        .default_value(60.0) // 60 seconds
        .nargs(1)
        .scan<'g', double>()
        .store_into(options.timelimit);

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
        .store_into(options.num_optimization_iters);

    // specific for potential methods
    program.add_argument("--potential_alpha")
        .help("how much weight to put on potential")
        .default_value(0.2)
        .nargs(1)
        .scan<'g', double>()
        .store_into(options.alpha);

    // specific for A1, A2 methods
    program.add_argument("--ax_nearest")
        .help("how relaxed to have filtering step in A1, A2 heuristics")
        .default_value(size_t{0})
        .nargs(1)
        .scan<'u', size_t>()
        .store_into(options.nearest);

    // only greedy potential ternary is implemented, so this is only valid with
    // single_shot and greedy_potential strategy, and
    program.add_argument("--ternary")
        .help("whether to optimize over {-1, 0, 1} rather than F_2")
        .default_value(false)
        .implicit_value(true)
        .store_into(options.ternary);
}

// for those arguments that cannot be filled in directly
void fill_options(argparse::ArgumentParser &program, slp::Options &options) {
    std::string search_method = program.get<std::string>("--search_method");
    std::string optimization_strategy =
        program.get<std::string>("--optimization_strategy");

    if (search_method == "greedy_potential") {
        options.strategy = slp::SearchStrategy::GreedyPotential;
    } else if (search_method == "backtrack_potential") {
        options.strategy = slp::SearchStrategy::BacktrackingPotential;
    } else if (search_method == "BP") {
        options.strategy = slp::SearchStrategy::BP;
    } else if (search_method == "paar1") {
        options.strategy = slp::SearchStrategy::Paar1;
    } else if (search_method == "RNBP") {
        options.strategy = slp::SearchStrategy::RNBP;
    } else if (search_method == "A1") {
        options.strategy = slp::SearchStrategy::A1;
    } else if (search_method == "A2") {
        options.strategy = slp::SearchStrategy::A2;
    } else if (search_method == "MIP") {
        options.strategy = slp::SearchStrategy::MIP;
    } else {
        throw std::invalid_argument(
            "received invalid argument for search method");
    }

    if (optimization_strategy == "framework") {
        options.optimization_strategy = slp::OptimizationStrategy::Framework;
    } else if (optimization_strategy == "single_shot") {
        options.optimization_strategy = slp::OptimizationStrategy::SingleShot;
    } else if (optimization_strategy == "repeat_random") {
        options.optimization_strategy = slp::OptimizationStrategy::RepeatRandom;
    } else {
        throw std::invalid_argument(
            "received invalid argument for optimization strategy");
    }
}

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    argparse::ArgumentParser program("SLP_CLI_solver_interface");
    slp::Options options;
    add_arguments(program, options);
    try {
        program.parse_args(argc, argv);
        fill_options(program, options);
        validate_options(options);
        if (options.skip_first) {
            int _;
            std::cin >> _;
        }

        // if {-1, 0, 1} matrix
        size_t m, n;
        if (options.ternary) {
            std::vector<std::vector<int>> matrix;
            read_ternary_matrix(m, n, matrix);
            int adds_before = get_ternary_naive_adds(matrix);
            auto [method, num_saved] =
                slp::ternary::gp::run_greedy_potential(matrix, m, options);
            int adds_after = adds_before - num_saved;

            std::cout << "-------------------------------------------\n";
            std::cout << "Naive additions: " << adds_before << "\n";
            std::cout << "After additions: " << adds_after << "\n";

            std::cout << "Method: \n";
            for (auto &[col1, col2, sigma] : method)
                std::cout << "col1: " << col1 << ", col2: " << col2
                          << ", sigma: " << sigma << "\n";

        } else {
            // default run
            std::vector<uint64_t> matrix;
            read_binary_matrix(m, n, matrix);
            slp::Z2Matrix G(matrix, m, n);

            slp::Result result = slp::gf2::run(G, options);
            std::cout << "-------------------------------------------\n";
            std::cout << "Naive additions: " << result.additions_before << "\n";
            std::cout << "After additions: " << result.additions_after << "\n";

            std::cout << "Method: \n";
            for (auto &[idx1, idx2] : result.method.additions)
                std::cout << idx1 << ", " << idx2 << "\n";

            std::cout << "Output indices: " << "\n";
            for (const size_t o : result.method.outputs)
                std::cout << o << "\n";
        }
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    return 0;
}
