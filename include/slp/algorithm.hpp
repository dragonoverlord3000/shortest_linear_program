#pragma once

#include "slp/types.hpp"

namespace slp::gf2 {
Result run(const Z2Matrix &G, const Options &options);
std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_greedy_potential(std::vector<uint64_t> &G, const Options &options);
std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_backtrack_potential(std::vector<uint64_t> &G, const Options &options);

std::vector<std::pair<std::size_t, std::size_t>>
run_boyar_peralta(const std::vector<uint64_t> &G, std::size_t m, std::size_t n,
                  const slp::Options &options);

AdditionMethod convert_potential_method(
    const std::vector<uint64_t> &G, std::size_t m, std::size_t n,
    std::vector<std::pair<std::size_t, std::size_t>> &potential_method);
} // namespace slp::gf2
