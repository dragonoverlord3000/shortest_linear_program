#include "slp/boyar_peralta/internal.hpp"

#include <iostream>

namespace slp::gf2::bp {

namespace {
    void step() {
    }
}

// returns method
// pair `p_i` in method means that element `i` is constructed by taking
// B[p_i[0]] xor B[p_i[1]], where B is the basis
std::vector<std::pair<size_t, size_t>> run_A1(const std::vector<uint64_t>& G, 
        size_t m, size_t n, const slp::Options &options) {
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
        if (options.verbose)
            std::cout << "BP Round #" << num_rounds << std::endl;

        step(basis, targets, dist, m, additions, s_targets_missing);
    }

    return additions;
}

} // namespace slp::gf2::bp
