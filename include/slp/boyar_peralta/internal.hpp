#pragma once

#include "slp/types.hpp"

namespace slp::gf2 {

std::vector<std::pair<size_t, size_t>>
run_boyar_peralta(const std::vector<uint64_t> &G, size_t m, size_t n,
                  const slp::Options &options);

AdditionMethod
convert_bp_method(std::vector<uint64_t> &G, size_t m, size_t n,
                  std::vector<std::pair<size_t, size_t>> &additions);
} // namespace slp::gf2
