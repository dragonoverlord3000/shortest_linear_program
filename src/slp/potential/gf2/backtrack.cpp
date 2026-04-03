#include "slp/potential/internal.hpp"

#include <algorithm>
#include <iostream>

namespace slp::gf2::gp {

// anonymous namespace -> _backtrack only usable here
namespace {
// returns number of savings possible for G using branch and bound
// note that the potential is an upper bound on the total additions savings we
// can have
void _backtrack(
    std::vector<uint64_t> &G, int current_saved, int potential,
    std::vector<std::pair<size_t, size_t>> &current_method,
    std::vector<std::pair<size_t, size_t>> &best_method,
    int &best_saved, const slp::Options &options, size_t lvl = 0) {
    if (current_saved + potential <= best_saved)
        return;

    if (current_saved > best_saved) {
        if (options.verbose) {
            std::cout << "New best saved: " << current_saved << std::endl;
            std::cout << "Potential(G): " << potential << std::endl;
            std::cout << "level at: " << lvl << std::endl;
        }
        best_saved = current_saved;
        best_method = current_method;
    }

    if (lvl == options.max_level || potential <= 0)
        return;

    // the tuples hold (-(num_saved + alpha * new_potential), new
    // current saved, new potential, col1, col2)
    std::vector<std::tuple<double, int, int, size_t, size_t>> moves;
    moves.reserve(G.size() * (G.size() - 1) / 2);
    for (size_t col1 = 0; col1 < G.size(); col1++) {
        for (size_t col2 = col1 + 1; col2 < G.size(); col2++) {
            uint64_t new_col = G[col1] & G[col2];
            if (std::popcount(new_col) < 2)
                continue;
            auto [num_saved, potential_diff] =
                evaluate_move(G, col1, col2, new_col);
            int new_potential = potential + potential_diff;
            int new_saved = current_saved + num_saved;
            if (new_saved + new_potential <= best_saved)
                continue;
            double neg_score = -(num_saved + options.alpha * new_potential);
            moves.push_back({neg_score, new_saved, new_potential, col1, col2});
        }
    }
    std::sort(moves.begin(), moves.end());

    for (auto &[_, new_current_saved, new_potential, col1, col2] : moves) {
        if (new_potential + new_current_saved <= best_saved)
            continue;
        current_method.push_back({col1, col2});
        apply_move(G, col1, col2);
        _backtrack(G, new_current_saved, new_potential, current_method,
                   best_method, best_saved, options, lvl + 1);
        undo_move(G, col1, col2);
        current_method.pop_back();
    }
}
} // namespace

// returns (method, num_saved)
std::pair<std::vector<std::pair<size_t, size_t>>, int>
run_backtrack_potential(std::vector<uint64_t> &G, const slp::Options &options) {
    int potential = get_potential(G);
    std::vector<std::pair<size_t, size_t>> current_method,
        best_method;
    int current_saved = 0, best_saved = 0;
    _backtrack(G, current_saved, potential, current_method, best_method,
               best_saved, options);
    return {best_method, best_saved};
}

} // namespace slp::gf2
