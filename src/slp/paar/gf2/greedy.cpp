// this corresponds to the Paar1 algorithm
#include "slp/types.hpp"

#include <stdint.h>
#include <vector>

namespace slp::gf2 {

std::vector<std::pair<std::size_t, std::size_t>>
run_paar1(std::vector<uint64_t> &G, const slp::Options &options) {
    std::vector<std::pair<std::size_t, std::size_t>> additions;

    std::size_t best_pc = 2;
    while (best_pc > 1) {
        best_pc = 1;
        uint64_t best_comb = 0;
        uint64_t best_i = 0;
        uint64_t best_j = 0;
        for (std::size_t i = 0; i < G.size(); i++) {
            for (std::size_t j = i + 1; j < G.size(); j++) {
                uint64_t comb = G[i] & G[j];
                std::size_t pc = std::popcount(comb);
                if (pc <= best_pc)
                    continue;

                best_pc = pc;
                best_comb = comb;
                best_i = i;
                best_j = j;
            }
        }

        if (best_pc <= 1)
            continue;
        additions.push_back({best_i, best_j});
        G[best_i] ^= best_comb;
        G[best_j] ^= best_comb;
        G.push_back(best_comb);
    }

    return additions;
}
} // namespace slp::gf2
