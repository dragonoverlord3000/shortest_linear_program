#include "slp/types.hpp"

namespace slp::gf2 {

std::pair<std::vector<Z2Matrix>, std::vector<PreprocStep>>
preprocess(const Z2Matrix &G, const Options &options);

Result
post_preprocess(const Z2Matrix &G, // the original G matrix before preprocessing
                const std::vector<Result> &results,
                std::vector<PreprocStep> &preproc_steps);
} // namespace slp::gf2
