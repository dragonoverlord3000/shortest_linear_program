#include "slp/algorithm.hpp"
#include "slp/potential/internal.hpp"
#include "slp/boyar_peralta/internal.hpp"
#include "slp/paar/internal.hpp"
#include "slp/types.hpp"

// For the modulo 2 algorithms
namespace slp::gf2 {

Result run(const Z2Matrix &_G, const Options &options) {
    size_t m = _G.m;
    size_t n = _G.n;
    std::vector<uint64_t> G = _G.matrix;

    Result result;

    // m is converted to int here, since we need its negation
    result.additions_before = gp::naive_additions(G, m);

    if (options.strategy == SearchStrategy::GreedyPotential) {
        auto [method, num_add_saved] = gp::run_greedy_potential(G, options);
        result.additions_after = result.additions_before - num_add_saved;
        result.method = gp::convert_potential_method(_G.matrix, m, n, method);
    } else if (options.strategy == SearchStrategy::BacktrackingPotential) {
        auto [method, num_add_saved] = gp::run_backtrack_potential(G, options);
        result.additions_after = result.additions_before - num_add_saved;
        result.method = gp::convert_potential_method(_G.matrix, m, n, method);
    } else if (options.strategy == SearchStrategy::BP) {
        std::vector<std::pair<size_t, size_t>> additions = bp::run_BP(G, m, n, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::RNBP) {
        std::vector<std::pair<size_t, size_t>> additions = bp::run_RNBP(G, m, n, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::A1) {
        std::vector<std::pair<size_t, size_t>> additions = bp::run_Ax(G, m, n, 1, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::A2) {
        std::vector<std::pair<size_t, size_t>> additions = bp::run_Ax(G, m, n, 2, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::Paar1) {
        std::vector<std::pair<size_t, size_t>> additions = paar::run_paar1(G, options);
        result.method = paar::convert_paar_method(_G.matrix, m, n, additions);
        result.additions_after = result.method.additions.size();
    }
    return result;
}

} // namespace slp::gf2
