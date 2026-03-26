#pragma once

#include <slp/types.hpp>

namespace slp::gf2 {

std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_greedy_boyar_peralta(std::vector<uint64_t> &G, std::size_t m, std::size_t n,
                         const slp::Options &options);

} // namespace slp::gf2
