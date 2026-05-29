#include "slp/boyar_peralta/internal.hpp"
#include "slp/utils/utils.hpp"

#include <chrono>
#include <iostream>
#include <random>

// this implementation is inspired by:
// https://github.com/thomaspeyrin/XORreduce/blob/master/main_globalopt.cpp

namespace slp::gf2::bp {
std::mt19937 rand_generator_rnbp;

namespace {
void step(Basis &basis, const std::vector<uint64_t> &targets,
          std::vector<size_t> &dist, size_t m,
          std::vector<std::pair<size_t, size_t>> &additions,
          std::unordered_set<uint64_t> &s_targets_missing,
          std::uniform_int_distribution<uint64_t> &rand_distribution) {
    size_t best_dist_norm = 0;
    size_t best_dist_sum = std::numeric_limits<size_t>::max();

    // DIST1 START
    // if dist[some_target] = 1 then make dist[some_target] = 0
    size_t dist1_idx = std::numeric_limits<size_t>::max();
    for (size_t i = 0; i < dist.size(); i++)
        if (dist[i] == 1) {
            dist1_idx = i;
            break;
        }

    if (dist1_idx != std::numeric_limits<size_t>::max()) {
        uint64_t new_b = targets[dist1_idx];
        for (size_t i = 0; i < basis.size(); i++) {
            if (!basis.contains(new_b ^ basis[i]))
                continue;
            for (size_t j = i + 1; j < basis.size(); j++)
                if (new_b == (basis[i] ^ basis[j])) {
                    std::vector<size_t> new_dist(m);
                    evaluate_move_bp(basis, targets, new_dist, dist, new_b);
                    apply_move_bp(basis, new_dist, dist, additions, i, j,
                                  new_b);
                    s_targets_missing.erase(new_b);
                    return;
                }
        }
        return;
    }
    // DIST1 END

    // the best possible distance values and corresponding candidate index pairs
    std::vector<std::vector<size_t>> best_dist;
    std::vector<std::pair<size_t, size_t>> candidates;

    for (size_t i = 0; i < basis.size(); i++) {
        for (size_t j = i + 1; j < basis.size(); j++) {
            uint64_t new_b = basis[i] ^ basis[j];
            if (basis.contains(new_b))
                continue;

            std::vector<size_t> new_dist(m);
            auto [cur_d, cur_nd] =
                evaluate_move_bp(basis, targets, new_dist, dist, new_b);
            if ((cur_d < best_dist_sum) ||
                (cur_d == best_dist_sum && cur_nd > best_dist_norm)) {
                best_dist_sum = cur_d;
                best_dist_norm = cur_nd;

                best_dist.clear();
                best_dist.push_back(new_dist);
                candidates.clear();
                candidates.push_back({i, j});
            } else if (cur_d == best_dist_sum && cur_nd == best_dist_norm) {
                candidates.push_back({i, j});
                best_dist.push_back(new_dist);
            }
        }
    }

    size_t rand_idx =
        rand_distribution(rand_generator_rnbp) % candidates.size();
    auto [best_i, best_j] = candidates[rand_idx];
    uint64_t best_b = basis[best_i] ^ basis[best_j];
    apply_move_bp(basis, best_dist[rand_idx], dist, additions, best_i, best_j,
                  best_b);
}
} // namespace

// returns method
// pair `p_i` in method means that element `i` is constructed by taking
// B[p_i[0]] xor B[p_i[1]], where B is the basis
std::vector<std::pair<size_t, size_t>> run_RNBP(const std::vector<uint64_t> &G,
                                                size_t m, size_t n,
                                                const slp::Options &options) {
    if (m == 0)
        return {};
    const auto deadline = std::chrono::steady_clock::now() +
                          std::chrono::duration<double>(options.timelimit);
    assert(m <= 64 && n <= 64);

    rand_generator_rnbp.seed(options.temp_seed ? options.temp_seed
                                               : options.seed);
    std::uniform_int_distribution<uint64_t> rand_distribution(
        0, std::numeric_limits<uint64_t>::max());

    // each row of G is a target, G[j] is a column (i.e. variable)
    Basis basis(n, options.reachable_strategy);
    std::unordered_set<uint64_t> s_targets_missing;
    std::vector<uint64_t> targets;
    std::vector<size_t> dist;
    init_bp(G, basis, s_targets_missing, targets, dist, m, n);

    std::vector<std::pair<size_t, size_t>> additions;
    int num_rounds = 0;
    while (!s_targets_missing.empty() && remaining_seconds(deadline) > 0.0) {
        num_rounds++;
        if (options.verbose) {
            std::cout << "RNBP Round #" << num_rounds << std::endl;
            std::cout << "Distance Vector: ";
            print_vec(dist);
        }

        step(basis, targets, dist, m, additions, s_targets_missing,
             rand_distribution);
    }

    // if time ran out before all could be computed, just do the rest naively
    fill_missing_additions_naive(s_targets_missing, n, additions);
    return additions;
}

} // namespace slp::gf2::bp
