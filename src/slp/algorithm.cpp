#include "slp/algorithm.hpp"
#include "slp/potential/internal.hpp"

// For the modulo 2 algorithms
namespace slp::gf2 {

Result run(const Z2Matrix &_G, const Options &options) {
    std::size_t m = _G.m;
    std::size_t n = _G.n;
    std::vector<uint64_t> G = _G.matrix;

    Result result;
    // m is converted to int here, since we need its negation
    result.additions_before = naive_additions(G, m);

    if (options.strategy == SearchStrategy::GreedyPotential) {
        auto [method, num_add_saved] = run_greedy_potential(G, options);
        result.additions_after = result.additions_before - num_add_saved;
        result.method = convert_potential_method(_G.matrix, m, n, method);
    } else if (options.strategy == SearchStrategy::BacktrackingPotential) {
        auto [method, num_add_saved] = run_backtrack_potential(G, options);
        result.additions_after = result.additions_before - num_add_saved;
        result.method = convert_potential_method(_G.matrix, m, n, method);
    }
    return result;
}

} // namespace slp::gf2
