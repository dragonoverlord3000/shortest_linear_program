#pragma once

#include "slp/types.hpp"

namespace slp::gf2 {

std::vector<std::pair<std::size_t, std::size_t>>
run_paar1(std::vector<uint64_t> &G, const slp::Options &options);

AdditionMethod convert_paar_method(
    const std::vector<uint64_t> &G, std::size_t m, std::size_t n,
    std::vector<std::pair<std::size_t, std::size_t>> &paar_method);
} // namespace slp::gf2
