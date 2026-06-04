#include <bits/stdc++.h>
using namespace std;

inline bool invertible9(array<uint16_t, 9> B) {
    int r = 0;
    for (int c = 0; c < 9; ++c) {
        int piv = -1;
        for (int i = r; i < 9; ++i) {
            if (B[i] & (1u << c)) {
                piv = i;
                break;
            }
        }
        if (piv == -1)
            continue;
        swap(B[r], B[piv]);
        for (int i = 0; i < 9; ++i) {
            if (i != r && (B[i] & (1u << c)))
                B[i] ^= B[r];
        }
        ++r;
    }
    return r == 9;
}

inline uint16_t rand9_biased(mt19937 &rng, double bit_p) {
    bit_p = clamp(bit_p, 0.0, 1.0);
    const uint32_t thresh = static_cast<uint32_t>(bit_p * 4294967296.0); // 2^32

    uint16_t x = 0;
    for (int b = 0; b < 9; ++b) {
        if (rng() < thresh)
            x |= static_cast<uint16_t>(1u << b);
    }
    return x; // 0..511
}

// smaller bit_p => sparser rows (heavier penalty for many 1-bits)
inline void fill_random_GL9(mt19937 &rng, array<uint16_t, 9> &B, double bit_p) {
    do {
        for (auto &row : B) {
            do {
                row = rand9_biased(rng, bit_p);
            } while (row == 0); // enforce [1, 511]
        }
    } while (!invertible9(B));
}
