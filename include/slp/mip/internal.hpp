#pragma once

#include "slp/types.hpp"

namespace slp::mip {

Result run_MIP(const Z2Matrix &G, const size_t m, const size_t n,
               const slp::Options &options);

} // namespace slp::mip
