#pragma once

#include "slp/types.hpp"

namespace slp::gf2 {
Result run(const Z2Matrix &G, const Options &options);
std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_greedy_potential(std::vector<uint64_t> &G, const Options &options);
std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_backtrack_potential(std::vector<uint64_t> &G, const Options &options);
} // namespace slp::gf2
