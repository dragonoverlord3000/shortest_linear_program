#include "slp/boyar_peralta/internal.hpp"
#include "slp/types.hpp"
#include "slp/utils/uitls.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

// this implementation is inspired by:
// https://github.com/rub-hgi/shorter_linear_slps_for_mds_matrices/blob/master/slp_heuristic.cpp
// https://eprint.iacr.org/2025/1493.pdf

namespace slp::gf2::bp {

namespace {

void step(Basis &basis, const std::vector<uint64_t> &targets,
          std::vector<size_t> &dist, size_t m,
          std::vector<std::pair<size_t, size_t>> &additions,
          std::unordered_set<uint64_t> &s_targets_missing) {
    size_t best_dist_norm = 0;
    size_t best_dist_sum = std::numeric_limits<size_t>::max();
    std::vector<size_t> best_dist(m);
    size_t best_i = 0;
    size_t best_j = 0;

    for (size_t i = 0; i < basis.size(); i++) {
        for (size_t j = i + 1; j < basis.size(); j++) {
            uint64_t new_b = basis[i] ^ basis[j];
            if (basis.contains(new_b))
                continue;
            // if dist[some_target] = 1 then make dist[some_target] = 0
            if (s_targets_missing.count(new_b)) {
                std::vector<size_t> new_dist(m);
                evaluate_move_bp(basis, targets, new_dist, dist, new_b);
                apply_move_bp(basis, new_dist, dist, additions, i, j, new_b);
                s_targets_missing.erase(new_b);
                return;
            }

            std::vector<size_t> new_dist(m);
            auto [cur_d, cur_nd] =
                evaluate_move_bp(basis, targets, new_dist, dist, new_b);
            if ((cur_d < best_dist_sum) ||
                (cur_d == best_dist_sum && cur_nd > best_dist_norm)) {
                best_dist_sum = cur_d;
                best_dist_norm = cur_nd;
                best_dist = new_dist;
                best_i = i;
                best_j = j;
            }
        }
    }

    uint64_t best_b = basis[best_i] ^ basis[best_j];
    apply_move_bp(basis, best_dist, dist, additions, best_i, best_j, best_b);
}
} // namespace

// returns method
// pair `p_i` in method means that element `i` is constructed by taking
// B[p_i[0]] xor B[p_i[1]], where B is the basis
std::vector<std::pair<size_t, size_t>> run_BP(const std::vector<uint64_t> &G,
                                              size_t m, size_t n,
                                              const slp::Options &options) {
    assert(m <= 64 && n <= 64);

    // each row of G is a target, G[j] is a column (i.e. variable)
    Basis basis(n, options.reachable_strategy);
    std::unordered_set<uint64_t> s_targets_missing;
    std::vector<uint64_t> targets;
    std::vector<size_t> dist;
    init_bp(G, basis, s_targets_missing, targets, dist, m, n);

    std::vector<std::pair<size_t, size_t>> additions;
    int num_rounds = 0;
    while (!s_targets_missing.empty()) {
        num_rounds++;
        if (options.verbose) {
            std::cout << "BP Round #" << num_rounds << std::endl;
            std::cout << "Distance Vector: ";
            print_vec(dist);
        }

        step(basis, targets, dist, m, additions, s_targets_missing);
    }

    return additions;
}
} // namespace slp::gf2::bp
