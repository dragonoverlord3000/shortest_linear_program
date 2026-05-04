#include "slp/potential/internal.hpp"

#include <algorithm>
#include <cassert>
#include <limits>

namespace {
std::pair<std::vector<int>, size_t>
get_new_col(const std::vector<std::vector<int>> &G, size_t col1, size_t col2,
            int sigma, size_t m) {
    std::vector<int> new_col(m);
    size_t cnt = 0;
    for (size_t i = 0; i < m; i++)
        if (G[i][col1] == sigma * G[i][col2] && G[i][col1] != 0) {
            new_col[i] = G[i][col1];
            cnt++;
        }
    return {new_col, cnt};
}

int saving_part(int t) { return std::max(0, t - 1); }

int contribution_vec_vec(const std::vector<int> &a, const std::vector<int> &b,
                         size_t m) {
    int same = 0;
    int opposite = 0;

    for (size_t i = 0; i < m; i++) {
        if (a[i] == 0 || b[i] == 0)
            continue;

        if (a[i] == b[i])
            same++;
        else if (a[i] == -b[i])
            opposite++;
    }

    return saving_part(same) + saving_part(opposite);
}

int contribution_vec_col(const std::vector<int> &a,
                         const std::vector<std::vector<int>> &G, size_t col,
                         size_t m) {
    int same = 0;
    int opposite = 0;

    for (size_t i = 0; i < m; i++) {
        int b = G[i][col];

        if (a[i] == 0 || b == 0)
            continue;

        if (a[i] == b)
            same++;
        else if (a[i] == -b)
            opposite++;
    }

    return saving_part(same) + saving_part(opposite);
}

int contribution_col_col(const std::vector<std::vector<int>> &G, size_t col1,
                         size_t col2, size_t m) {
    int same = 0;
    int opposite = 0;

    for (size_t i = 0; i < m; i++) {
        int a = G[i][col1];
        int b = G[i][col2];

        if (a == 0 || b == 0)
            continue;

        if (a == b)
            same++;
        else if (a == -b)
            opposite++;
    }

    return saving_part(same) + saving_part(opposite);
}

std::pair<int, int>
evaluate_move_ternary(const std::vector<std::vector<int>> &G, size_t col1,
                      size_t col2, size_t m, int sigma,
                      const std::vector<int> &new_col) {
    assert(col2 > col1);
    assert(!G.empty());

    size_t cols = G[0].size();

    std::vector<int> new_col1(m), new_col2(m);

    int support_size = 0;

    for (size_t i = 0; i < m; i++) {
        new_col1[i] = G[i][col1] - new_col[i];
        new_col2[i] = G[i][col2] - sigma * new_col[i];

        if (new_col[i] != 0)
            support_size++;
    }

    int saved = support_size - 1;

    int potential_diff = 0;

    // pair (col1, col2) changes
    potential_diff += contribution_vec_vec(new_col1, new_col2, m);
    potential_diff -= contribution_col_col(G, col1, col2, m);

    // new pairs (new column, col1) and (new column, col2)
    // these should usually be zero by construction, but computing them directly
    // makes the formula robust and keeps the evaluator definitionally correct
    potential_diff += contribution_vec_vec(new_col, new_col1, m);
    potential_diff += contribution_vec_vec(new_col, new_col2, m);

    for (size_t j = 0; j < cols; j++) {
        if (j == col1 || j == col2)
            continue;

        // old pairs involving col1 and col2 are replaced by new pairs
        potential_diff += contribution_vec_col(new_col1, G, j, m);
        potential_diff -= contribution_col_col(G, col1, j, m);

        potential_diff += contribution_vec_col(new_col2, G, j, m);
        potential_diff -= contribution_col_col(G, col2, j, m);

        // the newly appended column contributes with every old column
        potential_diff += contribution_vec_col(new_col, G, j, m);
    }

    return {saved, potential_diff};
}

void apply_move_ternary(std::vector<std::vector<int>> &G, size_t col1,
                        size_t col2, size_t m, int sigma) {
    // sanity check
    assert(col2 > col1);

    for (size_t i = 0; i < m; i++) {
        int val = 0;
        if (G[i][col1] == sigma * G[i][col2] && G[i][col1] != 0) {
            val = G[i][col1];
            G[i][col1] = 0;
            G[i][col2] = 0;
        }
        G[i].push_back(val);
    }
}
} // namespace

namespace slp::ternary::gp {
// returns (method, num_saved)
// the method includes a sign
std::pair<std::vector<std::tuple<size_t, size_t, int>>, int>
run_greedy_potential(std::vector<std::vector<int>> &G, size_t m,
                     const slp::Options &options) {
    assert(G.size() == m);
    assert(!G.empty());

    int total_saved = 0;
    std::vector<std::tuple<size_t, size_t, int>> method;

    // actual potential is not needed, since we just care about relative
    // ordering -> determined by the potential difference only int potential =
    // get_potential(G);

    while (true) {
        size_t cols = G.empty() ? 0 : G[0].size();
        int max_saved = 0;
        int best_saved = 0;
        // int best_potential_diff = 0;
        double best_score = std::numeric_limits<double>::lowest();
        size_t best_col1 = 0;
        size_t best_col2 = 0;
        int best_sigma = 0;

        for (size_t col1 = 0; col1 < cols; col1++) {
            for (size_t col2 = col1 + 1; col2 < cols; col2++) {
                for (int sigma : {-1, 1}) {
                    auto [new_col, support_size] =
                        get_new_col(G, col1, col2, sigma, m);
                    if (support_size < 2)
                        continue;

                    // simulate action
                    auto [saved, potential_diff] =
                        evaluate_move_ternary(G, col1, col2, m, sigma, new_col);
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
        apply_move_ternary(G, best_col1, best_col2, m, best_sigma);
        // potential += best_potential_diff;
        total_saved += best_saved;
        method.push_back({best_col1, best_col2, best_sigma});
    }

    return {method, total_saved};
}
} // namespace slp::gf2::ternary
