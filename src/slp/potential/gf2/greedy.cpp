#include "slp/potential/internal.hpp"
#include <limits>

namespace slp::gf2 {
// returns (method, num_saved)
std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_greedy_potential(std::vector<uint64_t> &G, const slp::Options& options) {
    int total_saved = 0;
    std::vector<std::pair<std::size_t, std::size_t>> method;

    int potential = get_potential(G);
    while (true) {
        int best_saved = 0;
        double best_score = std::numeric_limits<double>::lowest();
        std::size_t best_col1 = 0;
        std::size_t best_col2 = 0;

        for (std::size_t col1 = 0; col1 < G.size(); col1++) {
            for (std::size_t col2 = col1 + 1; col2 < G.size(); col2++) {
                int v = __builtin_popcount(G[col1] & G[col2]);
                if (v < 2)
                    continue;

                // simulate action
                auto [saved, potential_diff] = evaluate_move(G, col1, col2);
                best_saved = std::max(best_saved, saved);
                int potential_after = potential + potential_diff;
                double score = saved + options.alpha * potential_after;

                if (score > best_score) {
                    best_score = score;
                    best_col1 = col1;
                    best_col2 = col2;
                }
            }
        }
        // no more direct saving actions
        if (best_saved <= 0)
            break;
        // update G -> G'
        auto [saved, potential_diff] = evaluate_move(G, best_col1, best_col2);
        apply_move(G, best_col1, best_col2);
        potential += potential_diff;
        total_saved += saved;
        method.push_back({best_col1, best_col2});
    }

    return {method, total_saved};
}
} // namespace LAM::gf2
