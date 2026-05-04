#pragma once

#include "slp/types.hpp"

#include <cassert>
#include <stdint.h>
#include <vector>

namespace slp::gf2::gp {
int naive_additions(const std::vector<uint64_t> &G, int m);

int get_potential(const std::vector<uint64_t> &G);

std::pair<int, int> evaluate_move(const std::vector<uint64_t> &G, size_t col1,
                                  size_t col2, uint64_t new_col);

void apply_move(std::vector<uint64_t> &G, size_t col1, size_t col2);

void undo_move(std::vector<uint64_t> &G, size_t col1, size_t col2);

std::pair<std::vector<std::pair<size_t, size_t>>, int>
run_greedy_potential(std::vector<uint64_t> &G, const slp::Options &options);

std::pair<std::vector<std::pair<size_t, size_t>>, int>
run_backtrack_potential(std::vector<uint64_t> &G, const slp::Options &options);

AdditionMethod convert_potential_method(
    const std::vector<uint64_t> &G, size_t m, size_t n,
    std::vector<std::pair<size_t, size_t>> &potential_method);
} // namespace slp::gf2::gp

namespace slp::ternary::gp {
std::pair<std::vector<std::tuple<size_t, size_t, int>>, int>
run_greedy_potential(std::vector<std::vector<int>> &G, size_t m,
                     const slp::Options &options);
} // namespace slp::gf2::ternary
