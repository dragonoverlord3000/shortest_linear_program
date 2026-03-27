#include "slp/boyar_peralta/internal.hpp"
#include "slp/types.hpp"

#include <unordered_map>
#include <vector>

namespace slp::gf2 {
    // The `G` is the original `G` i.e. before applying `method`. Note that G[i] represents column i of G 
    // we are simply adding the target positions after the additions
    AdditionMethod convert_bp_method(std::vector<uint64_t>& G, std::size_t m, std::size_t n, std::vector<std::pair<std::size_t, std::size_t>>& additions) {
        assert(m <= 64);
        assert(n <= 64);
        AdditionMethod addition_method;
        addition_method.additions = additions;
        addition_method.outputs.assign(m, 0);

        std::vector<uint64_t> basis(n);
        for(std::size_t shift = 0; shift < n; shift++)
            basis[shift] = 1ULL << shift;

        std::unordered_map<uint64_t, std::vector<std::size_t>> targets2outputidx;
        for(std::size_t i = 0; i < m; i++) {
            uint64_t target = 0;
            for(std::size_t j = 0; j < n; j++) {
                if (G[j] & (1ULL << i)) 
                    target |= 1ULL << j;
            }
            if (!target) continue;
            targets2outputidx[target].push_back(i);
        }

        for(std::size_t i = 0; i < additions.size(); i++) {
            auto [l, r] = additions[i];
            uint64_t new_b = basis[l] ^ basis[r];
            if (targets2outputidx.count(new_b)) {
                for(std::size_t idx : targets2outputidx[new_b])
                    addition_method.outputs[idx] = basis.size();
                targets2outputidx.erase(new_b); // might not be necessary, should never create the same basis element twice
            }
            basis.push_back(new_b);
        }

        return addition_method;
    }
}
