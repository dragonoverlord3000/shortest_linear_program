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
    std::size_t m;
    std::size_t n;

    Z2Matrix(const std::vector<std::vector<int>> &_matrix) {
        assert(_matrix.size() != 0);
        m = _matrix.size();
        n = _matrix[0].size();

        assert(m <= 64);
        matrix.assign(n, 0);
        for (std::size_t i = 0; i < m; i++) {
            assert(_matrix[i].size() == n);
            for (std::size_t j = 0; j < n; j++) {
                assert(_matrix[i][j] == 0 || _matrix[i][j] == 1);
                matrix[j] |= _matrix[i][j] << i;
            }
        }
    }

    // Assumes each element of matrix represents a column
    Z2Matrix(std::vector<uint64_t>& matrix, std::size_t m, std::size_t n) : matrix(matrix), m(m), n(n) {
        assert(matrix.size() == n);
    };

    // e.g. for directly doing
    // >>> Z2Matrix G = {{1,1,1,0}, {0,1,1,1}}
    Z2Matrix(std::initializer_list<std::initializer_list<int>> init) {
        assert(init.size() != 0);

        m = init.size();
        n = init.begin()->size();

        assert(m <= 64);
        matrix.assign(n, 0);

        std::size_t i = 0;
        for (const std::initializer_list<int> &row : init) {
            assert(row.size() == n);

            std::size_t j = 0;
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

// TERNARY START
class TernaryMatrix {
  public:
    std::size_t m;
    std::size_t n;
    std::vector<std::vector<int>> matrix;

    TernaryMatrix(const std::vector<std::vector<int>> &_matrix) {
        assert(_matrix.size() != 0);
        m = _matrix.size();
        n = _matrix[0].size();

        matrix.assign(m, std::vector<int>(n));

        std::size_t i = 0;
        for (const std::vector<int> &row : _matrix) {
            assert(row.size() == n);

            std::size_t j = 0;
            for (const int x : row) {
                assert(std::abs(x) <= 1);
                matrix[i][j] = x;
                j++;
            }
            i++;
        }
    }
};
// TERNARY END

enum class SearchStrategy { GreedyPotential, BacktrackingPotential };

// Data passing helpers
struct Options {
    bool verbose = false;
    double alpha = 0.2;
    SearchStrategy strategy = SearchStrategy::GreedyPotential;
    std::size_t max_level = std::numeric_limits<std::size_t>::max();
};

struct Result {
    std::size_t additions_before = 0;
    std::size_t additions_after = 0;
    std::vector<std::pair<std::size_t, std::size_t>> transformation;
};
} // namespace slp
