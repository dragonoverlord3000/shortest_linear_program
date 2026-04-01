#pragma once

#include "slp/types.hpp"

namespace slp::gf2 {

std::vector<std::pair<size_t, size_t>>
run_paar1(std::vector<uint64_t> &G, const slp::Options &options);

AdditionMethod convert_paar_method(
    const std::vector<uint64_t> &G, size_t m, size_t n,
    std::vector<std::pair<size_t, size_t>> &paar_method);
} // namespace slp::gf2
