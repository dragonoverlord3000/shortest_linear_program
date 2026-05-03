#include "slp/potential/internal.hpp"

#include <limits>

namespace {
std::pair<std::vector<int>, size_t>
get_new_col(const std::vector<std::vector<int>> &G, size_t col1, size_t col2,
            int sigma, size_t m) {
    std::vector<int> new_col(m);
    size_t cnt = 0;
    for (size_t i = 0; i < m; i++)
        if (G[i][col1] == sigma * G[i][col2]) {
            new_col[i] = G[i][col1];
            cnt++;
        }
    return {new_col, cnt};
}

std::pair<int, int>
evaluate_move_ternary(const std::vector<std::vector<int>> &G, size_t col1,
                      size_t col2, int sigma, const std::vector<int> &new_col) {
    assert(col2 > col1);


}

void apply_move_ternary(std::vector<std::vector<int>> &G, size_t col1,
                        size_t col2) {}
} // namespace

namespace slp::gf2::ternary {
// returns (method, num_saved)
// the method includes a sign
std::pair<std::vector<std::tuple<size_t, size_t, int>>, int>
run_greedy_potential(std::vector<std::vector<int>> &G, size_t m, size_t n,
                     const slp::Options &options) {
    int total_saved = 0;
    std::vector<std::tuple<size_t, size_t, int>> method;

    // actual potential is not needed, since we just care about relative
    // ordering -> determined by the potential difference only int potential =
    // get_potential(G);

    while (true) {
        int max_saved = 0;
        int best_saved = 0;
        // int best_potential_diff = 0;
        double best_score = std::numeric_limits<double>::lowest();
        size_t best_col1 = 0;
        size_t best_col2 = 0;
        int best_sigma = 0;

        for (size_t col1 = 0; col1 < G.size(); col1++) {
            for (size_t col2 = col1 + 1; col2 < G.size(); col2++) {
                for (int sigma : {-1, 1}) {
                    auto [new_col, support_size] =
                        get_new_col(G, col1, col2, sigma, m);
                    if (support_size < 2)
                        continue;

                    // simulate action
                    auto [saved, potential_diff] =
                        evaluate_move_ternary(G, col1, col2, sigma, new_col);
                    if (saved <= 0)
                        continue; // we only consider saving moves

                    max_saved = std::max(max_saved, saved);

                    double score = saved + options.alpha * potential_diff;
                    if (score > best_score) {
                        best_score = score;

                        best_saved = saved;
                        // best_potential_diff = potential_diff;
                        best_col1 = col1;
                        best_col2 = col2;
                        best_sigma = sigma;
                    }
                }
            }
        }
        // no more direct saving actions
        if (max_saved <= 0)
            break;
        // update G -> G'
        apply_move_ternary(G, best_col1, best_col2);
        // potential += best_potential_diff;
        total_saved += best_saved;
        method.push_back({best_col1, best_col2, 1});
    }

    return {method, total_saved};
}
} // namespace slp::gf2::ternary
