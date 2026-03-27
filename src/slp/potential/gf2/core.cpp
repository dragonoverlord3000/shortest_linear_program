#include "slp/potential/internal.hpp"

#include <cassert>
#include <unordered_map>

namespace slp::gf2 {

// the additions used from simply carrying out the matrix multiplication with
// the standard O(n^3) algorithm G is m x n
int naive_additions(const std::vector<uint64_t> &G, int m) {
    int total = -m;
    for (const uint64_t &g : G)
        total += std::popcount(g);
    return total;
}

// get_potential directly -> the very first potential has to be computed
// directly, the rest are computed from the potential changes
int get_potential(const std::vector<uint64_t> &G) {
    int potential = 0;
    for (std::size_t col1 = 0; col1 < G.size(); col1++)
        for (std::size_t col2 = col1 + 1; col2 < G.size(); col2++) {
            int val = std::popcount(G[col1] & G[col2]);
            potential += val;
            if (val)
                potential--;
        }
    return potential;
}

// returns the integer pair: (# saved additions, the potential difference)
std::pair<int, int> evaluate_move(const std::vector<uint64_t> &G,
                                  std::size_t col1, std::size_t col2,
                                  uint64_t new_col) {
    // sanity check
    assert(col2 > col1);

    // The number of saved additions by combining col1 and col2 is the number of
    // additions we can combine i.e.
    // std::popcount(new_col), minus the cost of having to add new_col back
    // to col1 and col2 when computing
    int saved = -1 + std::popcount(new_col);

    uint64_t new_col1 = G[col1] ^ new_col;
    uint64_t new_col2 = G[col2] ^ new_col;

    // C(G,i,j) = -1 + std::popcount(new_col)
    int potential_diff = saved;

    // G[col1] -> new_col1, G[col2] -> new_col2
    for (std::size_t i = 0; i < G.size(); i++) {
        // add the new contribution of column 1 and column i, remove the old
        // contribution
        if (i != col1) {
            potential_diff +=
                std::popcount(G[i] & new_col1) - std::popcount(G[i] & G[col1]);
            potential_diff -= ((G[i] & new_col1) > 0) - ((G[i] & G[col1]) > 0);
        }
        // add the new contribution of column 2 and column i, remove the old
        // contribution
        if (i != col2) {
            potential_diff +=
                std::popcount(G[i] & new_col2) - std::popcount(G[i] & G[col2]);
            potential_diff -= ((G[i] & new_col2) > 0) - ((G[i] & G[col2]) > 0);
        }
        // note that C(G', col1, n+1) = C(G', col2, n+1) = C(G', col1, col2) = 0
        if (i != col1 && i != col2) {
            uint64_t t = G[i] & new_col;
            potential_diff += std::popcount(t);
            potential_diff -= (t > 0);
        }
    }

    return {saved, potential_diff};
}

void apply_move(std::vector<uint64_t> &G, int col1, int col2) {
    // sanity check
    assert(col2 > col1);

    uint64_t new_col = G[col1] & G[col2];
    uint64_t new_col1 = G[col1] ^ new_col;
    uint64_t new_col2 = G[col2] ^ new_col;

    // Update G -> G'
    G[col1] = new_col1;
    G[col2] = new_col2;
    G.push_back(new_col);
}

void undo_move(std::vector<uint64_t> &G, int col1, int col2) {
    // sanity check
    assert(col2 > col1);
    G[col1] |= G.back();
    G[col2] |= G.back();
    G.pop_back();
}

// CONVERT OUTPUT START
// The `G` is the original `G` i.e. before applying `method`. Note that G[i]
// represents column i of G This is essentially the converting the output to the
// convention of the Boyar-Peralta algorithm
AdditionMethod convert_potential_method(
    const std::vector<uint64_t> &G, std::size_t m, std::size_t n,
    std::vector<std::pair<std::size_t, std::size_t>> &potential_method) {
    // sanity check
    assert(G.size() == n);
    assert(m <= 64);
    assert(n <= 64);

    std::unordered_map<uint64_t, std::vector<std::size_t>> target2rows;
    std::vector<uint64_t> targets(m);
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < n; j++)
            if (G[j] & (1ULL << i))
                targets[i] |= 1ULL << j;
        target2rows[targets[i]].push_back(i);
    }

    AdditionMethod addition_method;
    addition_method.outputs.resize(m);

    // first form the basis elements explicitly required by potential method
    std::vector<uint64_t> basis;
    for (std::size_t i = 0; i < n; i++)
        basis.push_back(uint64_t{1} << i);
    for (std::pair<std::size_t, std::size_t> &p : potential_method) {
        addition_method.additions.push_back(p);
        uint64_t new_b = basis[p.first] ^ basis[p.second];
        basis.push_back(new_b);
        if (target2rows.count(new_b)) {
            for (std::size_t i : target2rows[new_b])
                addition_method.outputs[i] = basis.size() - 1;
            target2rows.erase(new_b);
        }
    }

    // form the remaining targets
    for (const std::pair<uint64_t, std::vector<std::size_t>> p : target2rows) {
        uint64_t target = p.first;
        std::size_t best_b_idx = 0;
        for (std::size_t idx = 1; idx < basis.size(); idx++)
            if (std::popcount(target ^ basis[idx]) <
                std::popcount(target ^ basis[best_b_idx]))
                best_b_idx = idx;

        std::size_t b_idx = best_b_idx;
        while (target ^ basis[b_idx]) {
            std::size_t tz = std::countr_zero(target ^ basis[b_idx]);
            uint64_t new_b = basis[tz] ^ basis[b_idx];
            basis.push_back(new_b);
            addition_method.additions.push_back({tz, b_idx});
            b_idx = basis.size() - 1;
        }

        for (std::size_t i : p.second)
            addition_method.outputs[i] = b_idx;
    }

    return addition_method;
}

// CONVERT OUTPUT END

} // namespace slp::gf2
