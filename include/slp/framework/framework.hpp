#include "slp/types.hpp"

#include <random>

namespace slp::gf2 {
std::tuple<Z2Matrix, std::vector<size_t>, std::vector<size_t>>
construct_new_G(const Z2Matrix &G, const Result &result, std::mt19937 &rng,
                const Options &options);

Result merge_results(const Z2Matrix &G, const Result &result,
                     const Z2Matrix &new_G, const Result &new_result,
                     const std::vector<size_t> &Si,
                     const std::vector<size_t> &So);
} // namespace slp::gf2
