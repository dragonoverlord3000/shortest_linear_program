#include "slp/algorithm.hpp"
#include "slp/types.hpp"

#include <chrono>
#include <iostream>
#include <random>
using namespace std;

namespace {
void print_G(const slp::Z2Matrix &G) {
    for (size_t i = 0; i < G.m; i++) {
        for (size_t j = 0; j < G.n; j++)
            cout << ((G.matrix[j] & (1ULL << i)) ? 1 : 0) << " ";
        cout << endl;
    }
}

void print_result(const slp::Result &result) {
    cout << "Additions before: " << result.additions_before << endl;
    cout << "Additions after: " << result.additions_after << endl;
    cout << "Transformation:" << endl;
    for (const std::pair<std::size_t, std::size_t> &t : result.method.additions)
        cout << "{" << t.first << ", " << t.second << "}" << endl;
    cout << "Output idxs: " << endl;
    for (const size_t o : result.method.outputs)
        cout << o << endl;
}
} // namespace

int main() {

    size_t n = 6;
    size_t num_tries = 2000;
    double prob = 0.6;

    slp::Options options;
    options.verbose = false;
    options.debug = true;
    options.strategy = slp::SearchStrategy::RNBP;
    options.alpha = 0.2;
    options.reachable_strategy = slp::ReachableStrategy::BacktracingSparseAware;
    options.nearest = 1;
    options.max_level = 1;
    options.optimization_strategy = slp::OptimizationStrategy::RepeatRandom;
    options.timelimit = 7.2;
    options.seed = 6283185;

    std::mt19937 rng;
    rng.seed(options.seed);
    std::bernoulli_distribution b_dist(prob);

    for (size_t iter = 0; iter < num_tries; iter++) {
        std::vector<std::vector<int>> matrix(2 * n + 1, vector<int>(2 * n));
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) {
                // upper left matrix (G1)
                matrix[i][j] = b_dist(rng);
                // lower right matrix (G2)
                matrix[i + n][j + n] = b_dist(rng);
            }
        // last row
        for (size_t i = 0; i < 2 * n; i++)
            matrix[2 * n][i] = b_dist(rng);

        cout << "--------------------------------------" << endl;
        cout << "Running with mixed row" << endl;

        slp::Z2Matrix G = matrix;
        std::cout << "G: " << std::endl;
        print_G(G);

        auto t0 = std::chrono::steady_clock::now();
        slp::Result result = slp::gf2::run(G, options);
        auto t1 = std::chrono::steady_clock::now();

        print_result(result);
        cout << "Duration: " << static_cast<std::chrono::nanoseconds>(t1 - t0)
             << endl;

        cout << "--------------------------------------" << endl;
        cout << "Running without mixed row" << endl;
        matrix.pop_back();
        G = matrix;
        std::cout << "G: " << std::endl;
        print_G(G);

        t0 = std::chrono::steady_clock::now();
        slp::Result result_wo = slp::gf2::run(G, options);
        t1 = std::chrono::steady_clock::now();

        print_result(result_wo);
        cout << "Duration: " << static_cast<std::chrono::nanoseconds>(t1 - t0)
             << endl;

        // sanity check
        cout << "Difference: "
             << result.additions_after - result_wo.additions_after << endl;
        if (result.additions_after < result_wo.additions_after)
            cout << "WITH IS BETTER: " << result.additions_after << ", "
                 << result_wo.additions_after << endl;
    }

    return 0;
}
