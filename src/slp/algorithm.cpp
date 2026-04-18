#include "slp/algorithm.hpp"
#include "slp/boyar_peralta/internal.hpp"
#include "slp/paar/internal.hpp"
#include "slp/potential/internal.hpp"
#include "slp/preprocess/preprocess.hpp"
#include "slp/types.hpp"

#include <iostream>

// for the modulo 2 algorithms
namespace slp::gf2 {

namespace {
Result run_framework(const Z2Matrix &_G, const Options &options) { return {}; }

Result run_heuristic(const Z2Matrix &_G, const Options &options) {
    size_t m = _G.m;
    size_t n = _G.n;
    std::vector<uint64_t> G(_G.matrix.begin(), _G.matrix.end());

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
        std::vector<std::pair<size_t, size_t>> additions =
            bp::run_BP(G, m, n, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::RNBP) {
        std::vector<std::pair<size_t, size_t>> additions =
            bp::run_RNBP(G, m, n, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::A1) {
        std::vector<std::pair<size_t, size_t>> additions =
            bp::run_Ax(G, m, n, 1, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::A2) {
        std::vector<std::pair<size_t, size_t>> additions =
            bp::run_Ax(G, m, n, 2, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::Paar1) {
        std::vector<std::pair<size_t, size_t>> additions =
            paar::run_paar1(G, options);
        result.method = paar::convert_paar_method(_G.matrix, m, n, additions);
        result.additions_after = result.method.additions.size();
    }

    return result;
}
} // namespace

Result run(const Z2Matrix &_G, const Options &options) {
    // TODO: add iterative refinement process (framework paper)
    std::vector<Z2Matrix> Gs;
    std::vector<PreprocStep> preproc_steps;

    if (options.use_preprocess) {
        std::pair<std::vector<Z2Matrix>, std::vector<PreprocStep>>
            preproc_G_step = preprocess(_G);
        Gs = preproc_G_step.first;
        preproc_steps = preproc_G_step.second;
    } else
        Gs = {_G};

    if (options.verbose && options.use_preprocess) {
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "preprocessing:" << std::endl;
        std::cout << "num matrix splits: " << Gs.size() << std::endl;
        std::cout << "num preproc steps: " << preproc_steps.size() << std::endl;
    }

    std::vector<Result> results;
    for (const Z2Matrix &G : Gs) {
        Result result;
        if (options.use_framework) {
            result = run_framework(G, options);
        } else {
            result = run_heuristic(G, options);
        }
        results.push_back(result);
    }

    Result result = post_preprocess(_G, results, preproc_steps);

    return result;
}

} // namespace slp::gf2
