#include "slp/types.hpp"

namespace slp::gf2 {

std::pair<std::vector<Z2Matrix>, std::vector<PreprocStep>>
preprocess(Z2Matrix &G);


Result post_preprocess(const std::vector<Z2Matrix> &G,
                       const std::vector<Result> &result,
                       std::vector<PreprocStep> &steps);

}
