#include "slp/preprocess/preprocess.hpp"
#include "slp/types.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <limits>

namespace slp::gf2 {

namespace {
constexpr size_t ZERO = std::numeric_limits<size_t>::max();
Z2Matrix row_preprocess(const Z2Matrix &G,
                        std::vector<PreprocStep> &preproc_steps) {
    size_t m = G.m;
    size_t n = G.n;
    std::vector<uint64_t> matrix = G.matrix;

    if (m == 0)
        return G;

    // row removal doesn't work if this is not the case
    assert(m <= 64 && m > 0);

    // triples (row, col, type) being the only set bit in row - or if type is
    // RowRemove0, then this is a zero row
    std::vector<std::tuple<size_t, size_t, PreprocType>> row_col_type;

    // find all row_col pairs we can remove
    for (size_t i = 0; i < m; i++) {
        size_t cnt = 0;
        size_t col = 0;
        for (size_t j = 0; j < n; j++) {
            if (!(matrix[j] & (1ULL << i)))
                continue;
            cnt++;
            col = j;
            if (cnt >= 2)
                break;
        }
        if (cnt == 0)
            row_col_type.push_back({i,
                                    ZERO, // used as temporary value
                                    PreprocType::RowRemove0});
        else if (cnt == 1)
            row_col_type.push_back({i, col, PreprocType::RowRemove1});
    }

    if (row_col_type.empty())
        return G;

    // remove rows -> going in reverse is easiest
    for (size_t i = 0; i < row_col_type.size(); i++) {
        auto [row, col, type] = row_col_type[row_col_type.size() - 1 - i];

        // record change
        PreprocStep preproc_step;
        preproc_step.type = type;
        preproc_step.rows.push_back(row);
        preproc_step.columns.push_back(col);
        preproc_steps.push_back(preproc_step);

        // remove row
        for (size_t j = 0; j < n; j++) {
            uint64_t low = matrix[j] & ((1ULL << row) - 1);
            uint64_t high = (matrix[j] >> 1) &
                            (((1ULL << (m - 1)) - 1) ^ ((1ULL << row) - 1));
            matrix[j] = low | high;
        }
    }

    Z2Matrix new_G(matrix, m - row_col_type.size(), n);
    return new_G;
}

Z2Matrix column_preprocess(const Z2Matrix &G,
                           std::vector<PreprocStep> &preproc_steps) {
    size_t m = G.m;
    size_t n = G.n;
    std::vector<uint64_t> matrix = G.matrix;

    // triples (row, column, PreprocType) for columns to be removed
    std::vector<std::tuple<size_t, size_t, PreprocType>> row_col_type;
    for (size_t j = 0; j < n; j++) {
        size_t pc = std::popcount(matrix[j]);
        if (pc == 0)
            row_col_type.push_back({ZERO, j, PreprocType::ColRemove0});
        else if (pc == 1)
            row_col_type.push_back(
                {std::countr_zero(matrix[j]), j, PreprocType::ColRemove1});
    }

    std::vector<PreprocStep> temp_preproc_steps;
    std::vector<uint64_t> new_matrix;
    size_t at = 0;
    for (size_t col = 0; col < n; col++) {
        if (at == row_col_type.size() || std::get<1>(row_col_type[at]) != col) {
            new_matrix.push_back(matrix[col]);
            continue;
        }
        auto [row, _, type] = row_col_type[at];
        at++;

        // record change
        PreprocStep preproc_step;
        preproc_step.type = type;
        preproc_step.rows.push_back(row);
        preproc_step.columns.push_back(col);
        temp_preproc_steps.push_back(preproc_step);
    }
    // add the column removals in reverse, because we are removing from the back
    std::reverse(temp_preproc_steps.begin(), temp_preproc_steps.end());
    for (const PreprocStep &preproc_step : temp_preproc_steps)
        preproc_steps.push_back(preproc_step);

    Z2Matrix new_G(new_matrix, m, n - row_col_type.size());
    return new_G;
}

/*
// note that separation does not produce any new rows and columns that are
// removable i.e. if there is a removable row or column in one of G1 or G2,
// then it was also there in G, so this is called last
std::pair<std::vector<Z2Matrix>, std::vector<PreprocStep>>
separate_preprocess() {
    return {{}, {}};
}
*/
} // namespace

std::pair<std::vector<Z2Matrix>, std::vector<PreprocStep>>
preprocess(const Z2Matrix &G) {

    std::vector<PreprocStep> preproc_steps;
    size_t prev_m = ZERO;
    size_t prev_n = ZERO;
    Z2Matrix G_prime = G;

    // removing columns can make some new rows removable
    while (prev_m != G_prime.m || prev_n != G_prime.n) {
        prev_m = G_prime.m;
        prev_n = G_prime.n;

        G_prime = row_preprocess(G_prime, preproc_steps);
        G_prime = column_preprocess(G_prime, preproc_steps);
    }

    return {{G_prime}, preproc_steps};
}

Result
post_preprocess(const Z2Matrix &G, // the original G matrix before preprocessing
                const std::vector<Result> &results,
                std::vector<PreprocStep> &preproc_steps) {
    Result result;
    result.additions_before = 0;
    result.additions_after = 0;

    // first combine separated matrices
    while (!preproc_steps.empty() &&
           preproc_steps.back().type == PreprocType::Separable) {
        assert(false); // not implemented yet
    }

    // this code is subject to change when separation preprocessing is added
    assert(results.size() == 1); // while separation not implemented
    result.additions_before += results[0].additions_before;
    result.additions_after += results[0].additions_after;
    result.method = results[0].method;

    // assumes remaining steps are adding back columns and rows
    size_t cur_n = G.n;
    for (const PreprocStep &preproc_step : preproc_steps)
        if (preproc_step.type == PreprocType::ColRemove0 ||
            preproc_step.type == PreprocType::ColRemove1)
            cur_n--;

    while (!preproc_steps.empty()) {
        PreprocStep preproc_step = preproc_steps.back();
        size_t row = preproc_step.rows.back();
        size_t col = preproc_step.columns.back();
        if (preproc_step.type == PreprocType::RowRemove0 ||
            preproc_step.type == PreprocType::RowRemove1) {
            // inserts before current row
            result.method.outputs.insert(result.method.outputs.begin() + row,
                                         col);
        } else if (preproc_step.type == PreprocType::ColRemove0) {
            cur_n++;
            // change the "name" (idx) of variables
            for (auto &[idx1, idx2] : result.method.additions) {
                if (idx1 >= col && idx1 != ZERO)
                    idx1++;
                if (idx2 >= col && idx2 != ZERO)
                    idx2++;
            }
            // change the outputs, since adding a column adds a dimension and
            // hence a new standard-basis to B
            for (size_t &o : result.method.outputs)
                if (o != ZERO && o >= col)
                    o++;
        } else if (preproc_step.type == PreprocType::ColRemove1) {
            cur_n++;
            // change "name" (idx) of variables
            for (auto &[idx1, idx2] : result.method.additions) {
                if (idx1 >= col && idx1 != ZERO)
                    idx1++;
                if (idx2 >= col && idx2 != ZERO)
                    idx2++;
            }
            // change the outputs, since adding a column adds a dimension and
            // hence a new standard-basis to B
            for (size_t &o : result.method.outputs)
                if (o != ZERO && o >= col)
                    o++;

            if (result.method.outputs[row] == ZERO) {
                result.method.outputs[row] = col;
            } else {
                // output[i] is the addition corresponding to
                result.method.additions.push_back(
                    {result.method.outputs[row], col});
                result.method.outputs[row] =
                    cur_n + result.method.additions.size() - 1;
                result.additions_after++;
            }
        } else
            assert(false); // no such type
        preproc_steps.pop_back();
    }

    return result;
}
} // namespace slp::gf2
