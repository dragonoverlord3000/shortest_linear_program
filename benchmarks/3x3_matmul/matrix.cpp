#include "matrix.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <stdint.h>

namespace matrix {
bool invertible9(std::array<uint16_t, 9> B) {
    int r = 0;
    for (int c = 0; c < 9; c++) {
        int piv = -1;
        for (int i = r; i < 9; i++) {
            if (B[i] & (1u << c)) {
                piv = i;
                break;
            }
        }
        if (piv == -1)
            continue;
        std::swap(B[r], B[piv]);
        for (int i = 0; i < 9; ++i) {
            if (i != r && (B[i] & (1u << c)))
                B[i] ^= B[r];
        }
        ++r;
    }
    return r == 9;
}

uint16_t rand9_biased(std::mt19937 &rng, double bit_p) {
    bit_p = std::clamp(bit_p, 0.0, 1.0);
    const uint32_t thresh = static_cast<uint32_t>(bit_p * 4294967296.0); // 2^32

    uint16_t x = 0;
    for (int b = 0; b < 9; b++) {
        if (rng() < thresh)
            x |= static_cast<uint16_t>(1u << b);
    }
    return x; // 0..511
}

// smaller bit_p => sparser rows (heavier penalty for many 1-bits)
void fill_random_GL9(std::mt19937 &rng, std::array<uint16_t, 9> &B,
                     double bit_p) {
    do {
        for (auto &row : B) {
            do {
                row = rand9_biased(rng, bit_p);
            } while (row == 0); // enforce [1, 511]
        }
    } while (!invertible9(B));
}

// BASIS CHANGE START
// G <- BG
inline void _change_basis_W(const std::vector<uint64_t> &G,
                            const std::array<uint16_t, 9> &B,
                            std::vector<uint64_t> &out) {
    out.assign(G.size(), 0);

    for (size_t k = 0; k < G.size(); k++) {
        uint64_t _g = G[k];
        uint64_t g = 0;

        // output bit i = parity(popcount(_g & B[i]))
        for (int i = 0; i < 9; i++) {
            unsigned x = static_cast<unsigned>(_g & B[i]);
            g |= static_cast<uint16_t>(__builtin_parity(x)) << i;
        }

        out[k] = g;
    }
}

// G <- GB
inline void _change_basis_UV(const std::vector<uint64_t> &G,
                             const std::array<uint16_t, 9> &B,
                             std::vector<uint64_t> &out) {
    out.assign(G.size(), 0);

    for (int col = 0; col < 9; col++)
        for (int k = 0; k < 9; k++)
            if (B[k] & (1 << col))
                out[col] ^= G[k];
}
// BASIS CHANGE END

} // namespace matrix
