#include "slp/types.hpp"

namespace slp::gf2::fw {

std::tuple<Z2Matrix, std::vector<size_t>, std::vector<size_t>>
construct_new_G(const Z2Matrix &G, const Result &result, const size_t gap,
                const size_t start, const Options& options);

Result merge_results(const Z2Matrix &G, const Result &result,
                     const Z2Matrix &new_G, const Result &new_result,
                     const std::vector<size_t> &Si,
                     const std::vector<size_t> &So);
} // namespace slp::gf2::fw
