#pragma once

#include <cassert>
#include <cstdlib>
#include <limits>
#include <stdint.h>
#include <vector>

namespace slp {

// MOD 2 START
class Z2Matrix {
  public:
    std::vector<uint64_t> matrix;
    size_t m;
    size_t n;

    Z2Matrix(const std::vector<std::vector<int>> &_matrix) {
        assert(_matrix.size() != 0);
        assert(_matrix[0].size() != 0);

        m = _matrix.size();
        n = _matrix[0].size();

        assert(m <= 64);
        assert(n <= 64);
        matrix.assign(n, 0);
        for (size_t i = 0; i < m; i++) {
            assert(_matrix[i].size() == n);
            for (size_t j = 0; j < n; j++) {
                assert(_matrix[i][j] == 0 || _matrix[i][j] == 1);
                matrix[j] |= _matrix[i][j] << i;
            }
        }
    }

    // Assumes each element of matrix represents a column
    Z2Matrix(const std::vector<uint64_t> &matrix, size_t m, size_t n)
        : matrix(matrix), m(m), n(n) {
        assert(matrix.size() == n);
        assert(n <= 64);
        assert(m <= 64);
    };

    // e.g. for directly doing
    // >>> Z2Matrix G = {{1,1,1,0}, {0,1,1,1}}
    Z2Matrix(std::initializer_list<std::initializer_list<int>> init) {
        assert(init.size() != 0);

        m = init.size();
        n = init.begin()->size();

        assert(m <= 64);
        assert(n <= 64);

        matrix.assign(n, 0);
        size_t i = 0;
        for (const std::initializer_list<int> &row : init) {
            assert(row.size() == n);

            size_t j = 0;
            for (int x : row) {
                assert(x == 0 || x == 1);
                matrix[j] |= static_cast<uint64_t>(x) << i;
                j++;
            }
            i++;
        }
    }
};
// MOD 2 END

enum class SearchStrategy {
    GreedyPotential,
    BacktrackingPotential,
    BP,
    RNBP,
    A1,
    A2,
    Paar1
};

// for boyar-peralta style algorithms
enum class ReachableStrategy {
    BruteForce,
    MITM,
    BacktracingSparseAware
}; // TODO: would like to add an adaptive option

enum class OptimizationStrategy { SingleShot, Framework, RepeatRandom };

// Data passing helpers
struct Options {
    bool verbose = false;
    // the heuristic to use
    SearchStrategy strategy = SearchStrategy::GreedyPotential;
    // Xiang based framework for optimization
    OptimizationStrategy optimization_strategy =
        OptimizationStrategy::SingleShot;
    size_t num_optimization_iters = std::numeric_limits<size_t>::max();
    double prob_framework_include =
        0.3; // TODO, find best parameter, or best scheduler
    double prob_random_restart = 0.05; // TODO, find best parameter

    // how many seconds to spend until stop
    double timelimit = 3600;

    // for potential based heuristics
    double alpha = 0.2;

    // for Boyar-Peralta based heuristics
    ReachableStrategy reachable_strategy =
        ReachableStrategy::BacktracingSparseAware;

    // for A1, A2, 0 is lowest, higher means more relaxed filtering step
    size_t nearest = 0;

    // for backtracking-based heuristics
    size_t max_level = std::numeric_limits<size_t>::max();

    // for preprocessing and postprocessing
    bool use_preprocess = true;
    bool use_postprocess = true;

    // for random methods
    uint64_t seed = 628318;
    uint64_t temp_seed = 0;
};

struct AdditionMethod {
    // we start with a basis B = [e_1, e_2, ..., e_n], each additions[k] = {i,j}
    // has 0 <= i < j < (k + n) and represents forming a new basis B_k via B_i
    // xor B_j
    std::vector<std::pair<size_t, size_t>> additions;
    // TODO: make dedicated struct for outputs entries (e.g. for ease of use,
    // when some rows are 0)
    std::vector<size_t> outputs; // outputs[i] is the index of the i'th output
                                 // vector in the basis formed by `additions`
};

struct Result {
    size_t additions_before = 0;
    size_t additions_after = 0;

    // we standardize all algorithms to have the same output as the Boyar
    // Peralta algorithm, as it is one of the simplest formats to read
    AdditionMethod method;
};

enum class PreprocType {
    RowRemove0,
    RowRemove1,
    ColRemove0,
    ColRemove1,
    Separable
};

struct PreprocStep {
    PreprocType type;
    std::vector<size_t> rows;
    std::vector<size_t> columns;

    // number of Gs that separate via separate preprocessing step
    size_t num_Gs = 0;
};
} // namespace slp
