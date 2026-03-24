#include <array>
#include <random>
#include <stdint.h>

namespace matrix {
bool invertible9(std::array<uint16_t, 9> B);
uint16_t rand9_biased(std::mt19937 &rng, double bit_p);
void fill_random_GL9(std::mt19937 &rng, std::array<uint16_t, 9> &B,
                     double bit_p);

void _change_basis_W(const std::vector<uint64_t> &G,
                            const std::array<uint16_t, 9> &B,
                            std::vector<uint64_t> &out);
void _change_basis_UV(const std::vector<uint64_t> &G,
                             const std::array<uint16_t, 9> &B,
                             std::vector<uint64_t> &out);
} // namespace matrix
