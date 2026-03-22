#include "slp/potential/internal.hpp"

namespace slp::gf2 {

// the additions used from simply carrying out the matrix multiplication with
// the standard O(n^3) algorithm G is m x n
int naive_additions(const std::vector<uint64_t> &G, int m) {
    int total = -m;
    for (const uint64_t &g : G)
        total += __builtin_popcount(g);
    return total;
}

// get_potential directly -> the very first potential has to be computed
// directly, the rest are computed from the potential changes
int get_potential(const std::vector<uint64_t> &G) {
    int potential = 0;
    for (std::size_t col1 = 0; col1 < G.size(); col1++)
        for (std::size_t col2 = col1 + 1; col2 < G.size(); col2++) {
            int val = __builtin_popcount(G[col1] & G[col2]);
            potential += val;
            if (val)
                potential--;
        }
    return potential;
}

// returns the integer pair: (# saved additions, the potential difference)
std::pair<int, int> evaluate_move(const std::vector<uint64_t> &G, std::size_t col1,
                                  std::size_t col2) {
    // sanity check
    assert(col2 > col1);

    uint64_t new_col = G[col1] & G[col2];

    // The number of saved additions by combining col1 and col2 is the number of
    // additions we can combine i.e.
    // __builtin_popcount(new_col), minus the cost of having to add new_col back
    // to col1 and col2 when computing
    int saved = -1 + __builtin_popcount(new_col);

    uint64_t new_col1 = G[col1] ^ new_col;
    uint64_t new_col2 = G[col2] ^ new_col;

    // C(G,i,j) = -1 + __builtin_popcount(new_col)
    int potential_diff = saved;

    // G[col1] -> new_col1, G[col2] -> new_col2
    for (std::size_t i = 0; i < G.size(); i++) {
        // add the new contribution of column 1 and column i, remove the old
        // contribution
        if (i != col1) {
            potential_diff += __builtin_popcount(G[i] & new_col1) -
                              __builtin_popcount(G[i] & G[col1]);
            potential_diff -= ((G[i] & new_col1) > 0) - ((G[i] & G[col1]) > 0);
        }
        // add the new contribution of column 2 and column i, remove the old
        // contribution
        if (i != col2) {
            potential_diff += __builtin_popcount(G[i] & new_col2) -
                              __builtin_popcount(G[i] & G[col2]);
            potential_diff -= ((G[i] & new_col2) > 0) - ((G[i] & G[col2]) > 0);
        }
        // note that C(G', col1, n+1) = C(G', col2, n+1) = C(G', col1, col2) = 0
        if (i != col1 && i != col2) {
            uint64_t t = G[i] & new_col;
            potential_diff += __builtin_popcount(t);
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

} // namespace LAM::gf2
