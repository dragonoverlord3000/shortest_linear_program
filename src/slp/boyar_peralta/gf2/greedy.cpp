#include <slp/boyar_peralta/internal.hpp>

#include <unordered_set>

// this implementation is inspired by:
// https://github.com/rub-hgi/shorter_linear_slps_for_mds_matrices/blob/master/slp_heuristic.cpp
// https://eprint.iacr.org/2025/1493.pdf

namespace slp::gf2 {
// returns (method, num_saved)
// pair `p_i` in method means that element `i` is constructed by taking B[p_i[0]] xor B[p_i[1]], where B is the basis
std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_greedy_boyar_peralta(std::vector<uint64_t> &G, std::size_t m, std::size_t n,
                         const slp::Options &options) {

    // each column of G is a target
    std::vector<uint64_t> basis;
    for(std::size_t shift = 0; shift < m; shift++)
        basis.push_back(1ll << shift);

    std::vector<std::size_t> dist(n);
    std::size_t targets_missing = n;
    for(std::size_t idx = 0; idx < n; idx++) {
        dist[idx] = std::popcount(G[idx]) - 1;
        if (!dist[idx]) targets_missing--;
    }

    int num_rounds = 0;
    std::vector<std::pair<std::size_t,std::size_t>> method;
    while (targets_missing) {
        num_rounds++;
        
    }

    std::size_t num_additions_required = 0;

    return {};
}
} // namespace slp::gf2
