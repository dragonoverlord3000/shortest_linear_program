#include "slp/preprocess/preprocess.hpp"
#include "slp/types.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>

namespace slp::gf2 {

namespace {
size_t get_naive_additions(const Z2Matrix &G) {
    size_t ans = 0;
    for (size_t i = 0; i < G.m; i++) {
        int cnt = -1;
        for (size_t j = 0; j < G.n; j++)
            if (G.matrix[j] & (1ULL << i))
                cnt++;
        ans += std::max(0, cnt);
    }

    return ans;
}

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

// note that separation does not produce any new rows and columns that are
// removable i.e. if there is a removable row or column in one of G1 or G2,
// then it was also there in G, so this is called last
std::vector<Z2Matrix>
separate_preprocess(const Z2Matrix &G,
                    std::vector<PreprocStep> &preproc_steps) {
    // setup nodes and edges in graph
    std::vector<std::pair<int, int>> nodes;
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> edges;
    for (size_t i = 0; i < G.m; i++) {
        bool prev_node = false;
        for (size_t j = 0; j < G.n; j++)
            if (G.matrix[j] & (1ULL << i)) {
                if (prev_node)
                    edges.push_back({{i, j}, nodes.back()});
                nodes.push_back({i, j});
                prev_node = true;
            }
    }
    if (nodes.empty())
        return {G};

    for (size_t j = 0; j < G.n; j++) {
        std::pair<int, int> prev_node = {-1, -1};
        for (size_t i = 0; i < G.m; i++)
            if (G.matrix[j] & (1ULL << i)) {
                if (prev_node.first != -1)
                    edges.push_back({prev_node, {i, j}});
                prev_node = {i, j};
            }
    }

    // find the connected components
    std::vector<Z2Matrix> Gs;
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> adj;
    for (auto &[node_1, node_2] : edges) {
        adj[node_1].push_back(node_2);
        adj[node_2].push_back(node_1);
    }
    std::set<std::pair<int, int>> seen;
    for (std::pair<int, int> &main_node : nodes) {
        if (seen.count(main_node))
            continue;
        std::queue<std::pair<int, int>> q;
        q.push(main_node);
        std::vector<std::pair<int, int>> new_G_set_bits;
        while (!q.empty()) {
            std::pair<int, int> node = q.front();
            q.pop();
            if (seen.count(node))
                continue;
            seen.insert(node);
            new_G_set_bits.push_back(node);
            if (!adj.count(node))
                continue;
            for (auto &neigh : adj[node])
                q.push(neigh);
        }

        // connected component to PreprocStep and matrix
        PreprocStep preproc_step;
        preproc_step.type = slp::PreprocType::Separable;
        std::set<int> unique_rows, unique_cols;
        for (auto &[r, c] : new_G_set_bits) {
            unique_rows.insert(r);
            unique_cols.insert(c);
        }
        preproc_step.rows.assign(unique_rows.begin(), unique_rows.end());
        preproc_step.columns.assign(unique_cols.begin(), unique_cols.end());

        std::unordered_map<int, int> oldrow2newrow, oldcol2newcol;
        int new_row = 0, new_col = 0;
        for (const int row : unique_rows)
            oldrow2newrow[row] = new_row++;
        for (const int col : unique_cols)
            oldcol2newcol[col] = new_col++;

        std::vector<uint64_t> matrix(new_col);
        for (auto &[row, col] : new_G_set_bits) {
            int r = oldrow2newrow.at(row);
            int c = oldcol2newcol.at(col);
            matrix[c] |= 1ULL << r;
        }
        Z2Matrix G(matrix, new_row, new_col);
        Gs.push_back(G);
        preproc_steps.push_back(preproc_step);
    }

    if (Gs.size() == 1)
        preproc_steps.pop_back();

    return Gs;
}
} // namespace

std::pair<std::vector<Z2Matrix>, std::vector<PreprocStep>>
preprocess(const Z2Matrix &G, const Options &options) {

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

    std::vector<Z2Matrix> Gs = {G_prime};
    if (options.use_separation)
        Gs = separate_preprocess(G_prime, preproc_steps);

    return {Gs, preproc_steps};
}

Result
post_preprocess(const Z2Matrix &G, // the original G matrix before preprocessing
                const std::vector<Result> &results,
                std::vector<PreprocStep> &preproc_steps) {
    Result result;
    result.additions_before = get_naive_additions(G);
    result.additions_after = 0;

    size_t sep_count = 0;
    size_t core_n = 0;
    for (int idx = (int)preproc_steps.size() - 1; idx >= 0; idx--) {
        if (preproc_steps[idx].type != PreprocType::Separable)
            break;
        core_n += preproc_steps[idx].columns.size();
        sep_count++;
    }

    bool has_separation = sep_count > 0;
    if (!has_separation) {
        assert(results.size() == 1);
        result.additions_after += results[0].additions_after;
        result.method = results[0].method;
    } else {
        assert(results.size() == sep_count);
        // first combine separated matrices
        size_t at = results.size() - 1;
        size_t offset = 0;
        std::vector<std::pair<size_t, size_t>> outputs;
        while (!preproc_steps.empty() &&
               preproc_steps.back().type == PreprocType::Separable) {
            PreprocStep preproc_step = preproc_steps.back();
            preproc_steps.pop_back();
            Result inner_result = results[at];

            // merge into result
            for (auto &[idx1, idx2] : inner_result.method.additions) {
                size_t i1 = idx1 < preproc_step.columns.size()
                                ? preproc_step.columns[idx1]
                                : idx1 + offset +
                                      (core_n - preproc_step.columns.size());
                size_t i2 = idx2 < preproc_step.columns.size()
                                ? preproc_step.columns[idx2]
                                : idx2 + offset +
                                      (core_n - preproc_step.columns.size());
                result.method.additions.push_back({i1, i2});
            }

            // sanity check
            assert(preproc_step.rows.size() ==
                   inner_result.method.outputs.size());

            for (size_t i = 0; i < inner_result.method.outputs.size(); i++) {
                size_t o = inner_result.method.outputs[i];
                size_t idx =
                    o < preproc_step.columns.size()
                        ? preproc_step.columns[o]
                        : o + offset + (core_n - preproc_step.columns.size());
                outputs.push_back({preproc_step.rows[i], idx});
            }

            offset += inner_result.method.additions.size();
            if (at)
                at--;

            result.additions_after += inner_result.additions_after;
        }
        std::sort(outputs.begin(), outputs.end());
        for (auto &[_, o] : outputs)
            result.method.outputs.push_back(o);
    }

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
