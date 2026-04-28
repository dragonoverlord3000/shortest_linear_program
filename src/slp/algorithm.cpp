#include "slp/algorithm.hpp"
#include "slp/boyar_peralta/internal.hpp"
#include "slp/framework/framework.hpp"
#include "slp/paar/internal.hpp"
#include "slp/postprocess/postprocess.hpp"
#include "slp/potential/internal.hpp"
#include "slp/preprocess/preprocess.hpp"
#include "slp/types.hpp"
#include "slp/utils/utils.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

// for the modulo 2 algorithms
namespace slp::gf2 {

namespace {

// little helper
template <class Duration>
double
remaining_seconds(const std::chrono::time_point<std::chrono::steady_clock,
                                                Duration> &deadline) {
    const auto now = std::chrono::steady_clock::now();
    return std::max(0.0, std::chrono::duration<double>(deadline - now).count());
}

// list of all search strategies to use in framework (note backtrack methods are
// too timeconsuming)
constexpr std::array<SearchStrategy, 6> all_search_strategies = {
    SearchStrategy::GreedyPotential,
    SearchStrategy::BP,
    SearchStrategy::RNBP,
    SearchStrategy::A1,
    SearchStrategy::A2,
    SearchStrategy::Paar1};

size_t get_naive_additions(const Z2Matrix &G, size_t m, size_t n) {
    size_t ans = 0;
    for (size_t i = 0; i < m; i++) {
        int cnt = -1;
        for (size_t j = 0; j < n; j++)
            cnt += (G.matrix[j] & (1ULL << i)) ? 1 : 0;
        ans += std::max(0, cnt);
    }
    return ans;
}

Result run_heuristic(const Z2Matrix &_G, const Options &options) {
    size_t m = _G.m;
    size_t n = _G.n;
    std::vector<uint64_t> G(_G.matrix.begin(), _G.matrix.end());

    Result result;

    // m is converted to int here, since we need its negation
    result.additions_before = get_naive_additions(_G, m, n);

    if (options.strategy == SearchStrategy::GreedyPotential) {
        if (options.verbose)
            std::cout << "Heuristic: Greedy Potential" << std::endl;
        auto [method, num_add_saved] = gp::run_greedy_potential(G, options);
        result.method = gp::convert_potential_method(_G.matrix, m, n, method);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::BacktrackingPotential) {
        if (options.verbose)
            std::cout << "Heuristic: Backtracking Potential" << std::endl;
        auto [method, num_add_saved] = gp::run_backtrack_potential(G, options);
        result.method = gp::convert_potential_method(_G.matrix, m, n, method);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::BP) {
        if (options.verbose)
            std::cout << "Heuristic: Boyar Peralta" << std::endl;
        std::vector<std::pair<size_t, size_t>> additions =
            bp::run_BP(G, m, n, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::RNBP) {
        if (options.verbose)
            std::cout << "Heuristic: Random Normal Boyar Peralta" << std::endl;
        std::vector<std::pair<size_t, size_t>> additions =
            bp::run_RNBP(G, m, n, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::A1) {
        if (options.verbose)
            std::cout << "Heuristic: A1" << std::endl;
        std::vector<std::pair<size_t, size_t>> additions =
            bp::run_Ax(G, m, n, 1, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::A2) {
        if (options.verbose)
            std::cout << "Heuristic: A2" << std::endl;
        std::vector<std::pair<size_t, size_t>> additions =
            bp::run_Ax(G, m, n, 2, options);
        result.method = bp::convert_bp_method(G, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else if (options.strategy == SearchStrategy::Paar1) {
        if (options.verbose)
            std::cout << "Heuristic: Paar1" << std::endl;
        std::vector<std::pair<size_t, size_t>> additions =
            paar::run_paar1(G, options);
        result.method = paar::convert_paar_method(_G.matrix, m, n, additions);
        result.additions_after = result.method.additions.size();
    } else {
        throw std::invalid_argument("Unsupported search strategy");
    }

    return result;
}

Result run_framework2(const Z2Matrix &_G, Options options,
                      const double timelimit) {
    const auto deadline = std::chrono::steady_clock::now() +
                          std::chrono::duration<double>(timelimit);

    std::mt19937 rng;
    rng.seed(options.seed);
    std::uniform_int_distribution<size_t> strategy_dist(
        0, all_search_strategies.size() - 1);
    std::uniform_int_distribution<uint64_t> temp_seed_dist(
        0, std::numeric_limits<uint64_t>::max());

    // get baseline
    options.temp_seed = temp_seed_dist(rng);
    options.strategy = all_search_strategies[strategy_dist(rng)];
    Result best_result = run_heuristic(_G, options);

    if (options.debug)
        validate_method(best_result.method, _G.n, "after run_heuristic 1");

    if (options.use_postprocess) {
        best_result.method = postprocess(best_result.method, _G.n);
        best_result.additions_after = best_result.method.additions.size();
    }

    if (options.debug)
        validate_method(best_result.method, _G.n, "after postprocess 1");

    // for small enough cases we can simply return them as is
    if (best_result.method.additions.size() < 3)
        return best_result;

    Result result = best_result;
    size_t iter = 0;
    while (remaining_seconds(deadline) > 0.0 &&
           iter < options.num_optimization_iters) {

        if (options.framework_restart) {
            // we do a restart here to get a new scheme to optimize
            options.temp_seed = temp_seed_dist(rng);
            options.strategy = all_search_strategies[strategy_dist(rng)];
            result = run_heuristic(_G, options);

            if (options.debug)
                validate_method(result.method, _G.n, "after run_heuristic 2");

            if (options.use_postprocess) {
                result.method = postprocess(result.method, _G.n);
                result.additions_after = result.method.additions.size();
            }
            if (options.debug)
                validate_method(result.method, _G.n, "after postprocess 2");
        }
        // for small enough cases we can simply return them as is
        if (result.method.additions.size() < 3)
            return result;

        size_t gap = result.method.additions.size() - 1;
        while (gap >= 2 && remaining_seconds(deadline) > 0.0) {
            bool improved_at_this_gap = false;

            size_t start = 0;
            while (gap >= 2 && start + gap <= result.method.additions.size() &&
                   remaining_seconds(deadline) > 0.0) {
                if (iter >= options.num_optimization_iters)
                    break;
                iter++;

                options.temp_seed = temp_seed_dist(rng);
                options.strategy = all_search_strategies[strategy_dist(rng)];

                // create the new G
                auto [new_G, Si, So] =
                    fw::construct_new_G(_G, result, gap, start, options);

                std::vector<Z2Matrix> Gs;
                std::vector<PreprocStep> preproc_steps;

                if (options.use_preprocess) {
                    std::pair<std::vector<Z2Matrix>, std::vector<PreprocStep>>
                        preproc_G_step = preprocess(new_G);
                    Gs = preproc_G_step.first;
                    preproc_steps = preproc_G_step.second;
                } else
                    Gs = {new_G};

                // optimize the new G
                std::vector<Result> results;
                for (const Z2Matrix &G : Gs) {
                    Result t_result;
                    t_result = run_heuristic(G, options);

                    if (options.debug)
                        validate_method(t_result.method, G.n,
                                        "after run_heuristic 3");

                    if (options.use_postprocess) {
                        t_result.method = postprocess(t_result.method, G.n);
                        t_result.additions_after =
                            t_result.method.additions.size();
                    }

                    if (options.debug)
                        validate_method(t_result.method, G.n,
                                        "after postprocess 3");
                    results.push_back(t_result);
                }

                Result new_result =
                    post_preprocess(new_G, results, preproc_steps);

                if (options.debug)
                    validate_method(new_result.method, new_G.n,
                                    "after post_preprocess 1");

                // merge new_G optimization back into G
                Result candidate =
                    fw::merge_results(_G, result, new_G, new_result, Si, So);

                if (options.debug)
                    validate_method(candidate.method, _G.n, "after merge 1");

                if (options.use_postprocess)
                    candidate.method = postprocess(candidate.method, _G.n);
                candidate.additions_after = candidate.method.additions.size();

                if (options.debug)
                    validate_method(candidate.method, _G.n,
                                    "after merge postprocess 1");

                if (candidate.method.additions.size() <
                    result.method.additions.size()) {
                    result = candidate;

                    if (result.method.additions.size() <
                        best_result.method.additions.size()) {
                        best_result = result;
                    }
                    if (result.method.additions.size() < 3)
                        break;

                    gap = result.method.additions.size() - 1;
                    start = 0;
                    improved_at_this_gap = true;
                } else
                    start++;
            }

            if (!improved_at_this_gap)
                gap--;
        }

        if (result.method.additions.size() <
            best_result.method.additions.size())
            best_result = result;
    }

    return best_result;
}

/*
Result run_framework(const Z2Matrix &_G, Options options,
                     const double timelimit) {
    const auto deadline = std::chrono::steady_clock::now() +
                          std::chrono::duration<double>(timelimit);

    std::mt19937 rng;
    rng.seed(options.seed);
    std::uniform_real_distribution<double> random_restart_dist(0, 1);
    std::uniform_int_distribution<size_t> strategy_dist(
        0, all_search_strategies.size() - 1);
    std::uniform_int_distribution<uint64_t> temp_seed_dist(
        0, std::numeric_limits<uint64_t>::max());

    // start by getting the base result
    options.temp_seed = temp_seed_dist(rng);
    options.strategy = all_search_strategies[strategy_dist(rng)];
    Result result = run_heuristic(_G, options);
    if (options.use_postprocess) {
        result.method = postprocess(result.method, _G.n);
        result.additions_after = result.method.additions.size();
    }

    Result best_result = result;

    for (size_t iter = 1; iter < options.num_optimization_iters &&
                          remaining_seconds(deadline) > 0.0;
         iter++) {
        options.temp_seed = temp_seed_dist(rng);
        options.strategy = all_search_strategies[strategy_dist(rng)];

        // random restart
        if (random_restart_dist(rng) < options.prob_random_restart) {
            result = run_heuristic(_G, options);
            if (options.use_postprocess) {
                result.method = postprocess(result.method, _G.n);
                result.additions_after = result.method.additions.size();
            }
            if (result.additions_after < best_result.additions_after)
                best_result = result;
        }

        // create the new G
        auto [new_G, Si, So] = fw::construct_new_G(_G, result, rng, options);

        std::vector<Z2Matrix> Gs;
        std::vector<PreprocStep> preproc_steps;

        if (options.use_preprocess) {
            std::pair<std::vector<Z2Matrix>, std::vector<PreprocStep>>
                preproc_G_step = preprocess(new_G);
            Gs = preproc_G_step.first;
            preproc_steps = preproc_G_step.second;
        } else
            Gs = {new_G};

        // optimize the new G
        std::vector<Result> results;
        for (const Z2Matrix &G : Gs) {
            Result t_result;
            t_result = run_heuristic(G, options);
            if (options.use_postprocess) {
                t_result.method = postprocess(t_result.method, G.n);
                t_result.additions_after = t_result.method.additions.size();
            }
            results.push_back(t_result);
        }

        Result new_result = post_preprocess(new_G, results, preproc_steps);

        // merge new_G optimization back into G
        result = fw::merge_results(_G, result, new_G, new_result, Si, So);
        if (options.use_postprocess)
            result.method = postprocess(result.method, _G.n);
        result.additions_after = result.method.additions.size();

        if (result.additions_after < best_result.additions_after)
            best_result = result;
    }

    return best_result;
}
*/

Result run_repeat_random(const Z2Matrix &_G, Options options,
                         const double timelimit) {
    const auto deadline = std::chrono::steady_clock::now() +
                          std::chrono::duration<double>(timelimit);

    std::mt19937 rng;
    rng.seed(options.seed);
    std::uniform_int_distribution<uint64_t> temp_seed_dist(
        0, std::numeric_limits<uint64_t>::max());

    // start by getting the default result
    Result result = run_heuristic(_G, options);
    if (options.use_postprocess) {
        result.method = postprocess(result.method, _G.n);
        result.additions_after = result.method.additions.size();
    }

    Result best_result = result;

    for (size_t iter = 1; iter < options.num_optimization_iters &&
                          remaining_seconds(deadline) > 0.0;
         iter++) {

        options.temp_seed = temp_seed_dist(rng);

        result = run_heuristic(_G, options);
        if (options.use_postprocess) {
            result.method = postprocess(result.method, _G.n);
            result.additions_after = result.method.additions.size();
        }
        if (result.additions_after < best_result.additions_after)
            best_result = result;
    }

    return best_result;
}

} // namespace

Result run(const Z2Matrix &_G, const Options &options) {
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
        // note that we split the time, but not max num iters, as the
        // sub-problems should be that much faster to complete
        double amt_time = options.timelimit / Gs.size();

        Result result;
        if (options.optimization_strategy ==
            slp::OptimizationStrategy::Framework) {
            result = run_framework2(G, options, amt_time);
        } else if (options.optimization_strategy ==
                   slp::OptimizationStrategy::SingleShot) {
            result = run_heuristic(G, options);
            if (options.use_postprocess) {
                result.method = postprocess(result.method, G.n);
                result.additions_after = result.method.additions.size();
            }
        } else if (options.optimization_strategy ==
                   slp::OptimizationStrategy::RepeatRandom) {
            result = run_repeat_random(G, options, amt_time);
        } else {
            throw std::invalid_argument(
                "Unsupported optimization strategy strategy");
        }
        results.push_back(result);
    }

    Result result = post_preprocess(_G, results, preproc_steps);

    return result;
}

} // namespace slp::gf2
