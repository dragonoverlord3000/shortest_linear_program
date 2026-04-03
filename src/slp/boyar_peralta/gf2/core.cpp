#include "slp/boyar_peralta/internal.hpp"
#include "slp/types.hpp"

#include <unordered_map>
#include <vector>

namespace slp::gf2::bp {

// initialize elements required for bp-type algorithms
void init_bp(const std::vector<uint64_t> &G, Basis &basis,
             std::unordered_set<uint64_t> &s_targets_missing,
             std::vector<uint64_t> &targets, std::vector<size_t> &dist,
             size_t m, size_t n) {

    for (size_t shift = 0; shift < n; shift++) {
        uint64_t b = 1ULL << shift;
        basis.add_element(b);
    }

    targets.assign(m, 0);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            if (G[j] & (1ULL << i))
                targets[i] |= 1ULL << j;
        }
        if (targets[i] == 0 || basis.contains(targets[i]))
            continue;
        s_targets_missing.insert(targets[i]);
    }

    dist.assign(m, 0);
    for (size_t idx = 0; idx < m; idx++) {
        if (targets[idx] == 0 || basis.contains(targets[idx])) {
            dist[idx] = 0;
            continue;
        }
        size_t d = 0;
        for (size_t col = 0; col < n; col++)
            d += (G[col] & (1ULL << idx)) > 0;
        dist[idx] = d - 1;
    }
}

// The `G` is the original `G` i.e. before applying `method`. Note that G[i]
// represents column i of G we are simply adding the target positions after the
// additions
AdditionMethod
convert_bp_method(std::vector<uint64_t> &G, size_t m, size_t n,
                  std::vector<std::pair<size_t, size_t>> &additions) {
    assert(m <= 64);
    assert(n <= 64);
    AdditionMethod addition_method;
    addition_method.additions = additions;
    addition_method.outputs.assign(m, 0);

    std::vector<uint64_t> basis(n);
    for (size_t shift = 0; shift < n; shift++)
        basis[shift] = 1ULL << shift;

    std::unordered_map<uint64_t, std::vector<size_t>> targets2outputidx;
    for (size_t i = 0; i < m; i++) {
        uint64_t target = 0;
        for (size_t j = 0; j < n; j++) {
            if (G[j] & (1ULL << i))
                target |= 1ULL << j;
        }
        if (!target)
            continue;
        targets2outputidx[target].push_back(i);
    }

    for (size_t i = 0; i < additions.size(); i++) {
        auto [l, r] = additions[i];
        uint64_t new_b = basis[l] ^ basis[r];
        if (targets2outputidx.count(new_b)) {
            for (size_t idx : targets2outputidx[new_b])
                addition_method.outputs[idx] = basis.size();
            targets2outputidx.erase(
                new_b); // might not be necessary, should never create the same
                        // basis element twice
        }
        basis.push_back(new_b);
    }

    return addition_method;
}
} // namespace slp::gf2::bp
