#include "slp/algorithm.hpp"
#include "slp/types.hpp"

#include <chrono>
#include <iostream>
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
    cout << "-----------------------------------------" << endl;
    cout << "Greedy Potential" << endl;
    slp::Z2Matrix G = {{1,1,1,0}, {0,1,1,1}};

    print_G(G);

    slp::Options options;
    options.verbose = false;
    options.strategy = slp::SearchStrategy::MIP;
    options.optimization_strategy = slp::OptimizationStrategy::SingleShot;
    options.timelimit = 1200;
    options.seed = 628;

    auto t0 = std::chrono::steady_clock::now();
    slp::Result result = slp::gf2::run(G, options);
    auto t1 = std::chrono::steady_clock::now();

    print_result(result);
    cout << "Duration: " << static_cast<std::chrono::nanoseconds>(t1 - t0)
         << endl;

    return 0;
}
