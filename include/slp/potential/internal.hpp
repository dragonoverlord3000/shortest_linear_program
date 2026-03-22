#pragma once

#include "slp/types.hpp"
#include <cassert>
#include <stdint.h>
#include <vector>

namespace slp::gf2 {
int naive_additions(const std::vector<uint64_t> &G, int m);

int get_potential(const std::vector<uint64_t> &G);

std::pair<int, int> evaluate_move(const std::vector<uint64_t> &G,
                                  std::size_t col1, std::size_t col2);

void apply_move(std::vector<uint64_t> &G, int col1, int col2);

void undo_move(std::vector<uint64_t> &G, int col1, int col2);

std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_greedy_potential(std::vector<uint64_t> &G, const slp::Options& options);

std::pair<std::vector<std::pair<std::size_t, std::size_t>>, int>
run_backtrack_potential(std::vector<uint64_t> &G, const slp::Options& options);
} // namespace slp::gf2
