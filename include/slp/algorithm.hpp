#pragma once

#include "slp/types.hpp"

namespace slp::gf2 {
Result run_heuristic(const Z2Matrix &G, const Options &options);

Result run(const Z2Matrix &G, const Options &options);
} // namespace slp::gf2
