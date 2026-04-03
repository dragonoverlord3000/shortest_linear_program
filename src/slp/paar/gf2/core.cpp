#include "slp/paar/internal.hpp"

#include <cassert>
#include <unordered_map>

namespace slp::gf2::paar {

    // basically the same as the potential method conversion
AdditionMethod convert_paar_method(
    const std::vector<uint64_t> &G, size_t m, size_t n,
    std::vector<std::pair<size_t, size_t>> &paar_method) {
    // sanity check
    assert(G.size() == n);
    assert(m <= 64);
    assert(n <= 64);

    std::unordered_map<uint64_t, std::vector<size_t>> target2rows;
    std::vector<uint64_t> targets(m);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++)
            if (G[j] & (1ULL << i))
                targets[i] |= 1ULL << j;
        target2rows[targets[i]].push_back(i);
    }

    AdditionMethod addition_method;
    addition_method.outputs.resize(m);

    // first form the basis elements explicitly required by potential method
    std::vector<uint64_t> basis;
    for (size_t i = 0; i < n; i++)
        basis.push_back(uint64_t{1} << i);
    for (std::pair<size_t, size_t> &p : paar_method) {
        addition_method.additions.push_back(p);
        uint64_t new_b = basis[p.first] ^ basis[p.second];
        basis.push_back(new_b);
        if (target2rows.count(new_b)) {
            for (size_t i : target2rows[new_b])
                addition_method.outputs[i] = basis.size() - 1;
            target2rows.erase(new_b);
        }
    }

    // form the remaining targets
    for (const std::pair<uint64_t, std::vector<size_t>> p : target2rows) {
        uint64_t target = p.first;
        size_t best_b_idx = 0;
        for (size_t idx = 1; idx < basis.size(); idx++)
            if (std::popcount(target ^ basis[idx]) <
                std::popcount(target ^ basis[best_b_idx]))
                best_b_idx = idx;

        size_t b_idx = best_b_idx;
        while (target ^ basis[b_idx]) {
            size_t tz = std::countr_zero(target ^ basis[b_idx]);
            uint64_t new_b = basis[tz] ^ basis[b_idx];
            basis.push_back(new_b);
            addition_method.additions.push_back({tz, b_idx});
            b_idx = basis.size() - 1;
        }

        for (size_t i : p.second)
            addition_method.outputs[i] = b_idx;
    }

    return addition_method;
}
} // namespace slp::gf2
