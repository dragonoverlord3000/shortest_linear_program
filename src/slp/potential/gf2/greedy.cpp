#include "slp/potential/internal.hpp"

#include <limits>

namespace slp::gf2::gp {
// returns (method, num_saved)
std::pair<std::vector<std::pair<size_t, size_t>>, int>
run_greedy_potential(std::vector<uint64_t> &G, const slp::Options &options) {
    int total_saved = 0;
    std::vector<std::pair<size_t, size_t>> method;

    // actual potential is not needed, since we just care about relative ordering -> determined by the potential difference only
    // int potential = get_potential(G);

    while (true) {
        int max_saved = 0;
        int best_saved = 0;
        // int best_potential_diff = 0;
        double best_score = std::numeric_limits<double>::lowest();
        size_t best_col1 = 0;
        size_t best_col2 = 0;

        for (size_t col1 = 0; col1 < G.size(); col1++) {
            for (size_t col2 = col1 + 1; col2 < G.size(); col2++) {
                uint64_t new_col = G[col1] & G[col2];
                int v = std::popcount(new_col);
                if (v < 2)
                    continue;

                // simulate action
                auto [saved, potential_diff] =
                    evaluate_move(G, col1, col2, new_col);
                if (saved <= 0) continue; // we only consider saving moves

                max_saved = std::max(max_saved, saved);
                // note that potential is constant within this loop, so for positive alpha 
                // alpha * (potential + potential_diff) has same ordering as alpha * potential_diff
                double score = saved + options.alpha * potential_diff; 

                if (score > best_score) {
                    best_score = score;

                    best_saved = saved;
                    // best_potential_diff = potential_diff;
                    best_col1 = col1;
                    best_col2 = col2;
                }
            }
        }
        // no more direct saving actions
        if (max_saved <= 0)
            break;
        // update G -> G'
        apply_move(G, best_col1, best_col2);
        // potential += best_potential_diff;
        total_saved += best_saved;
        method.push_back({best_col1, best_col2});
    }

    return {method, total_saved};
}
} // namespace slp::gf2
