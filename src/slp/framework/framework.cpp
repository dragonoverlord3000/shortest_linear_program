#include "slp/framework/framework.hpp"
#include "slp/types.hpp"
#include "slp/utils/utils.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>

namespace slp::gf2 {
std::tuple<Z2Matrix, std::vector<size_t>, std::vector<size_t>>
construct_new_G(const Z2Matrix &G, const Result &result, std::mt19937 &rng,
                const Options &options) {
    assert(G.n <= 64 && G.m <= 64);
    std::bernoulli_distribution bernoulli_dist(options.prob_framework_include);

    std::vector<uint64_t> basis(G.n + result.method.additions.size());
    for (uint64_t shift = 0; shift < G.n; shift++)
        basis[shift] = 1ULL << shift;

    size_t N = G.n + result.method.additions.size();
    std::vector<std::vector<size_t>> adj(N), r_adj(N);
    for (size_t i = 0; i < result.method.additions.size(); i++) {
        size_t idx = i + G.n;
        auto &[idx1, idx2] = result.method.additions[i];
        basis[idx] = basis[idx1] ^ basis[idx2];
        adj[idx1].push_back(idx);
        adj[idx2].push_back(idx);
        r_adj[idx].push_back(idx1);
        r_adj[idx].push_back(idx2);
    }
    std::unordered_set<size_t> output_idxs;
    for (const size_t &o : result.method.outputs)
        output_idxs.insert(o);

    // construct the unprocessed set So
    std::unordered_set<size_t> So_set;
    for (size_t idx = G.n; idx < G.n + result.method.additions.size(); idx++)
        if (bernoulli_dist(rng))
            So_set.insert(idx);

    // all the nodes used to form So
    std::unordered_set<size_t> So_group;
    std::queue<size_t> q;
    for (const size_t &idx : So_set)
        q.push(idx);
    while (!q.empty()) {
        size_t idx = q.front();
        q.pop();
        if (So_group.count(idx))
            continue;
        So_group.insert(idx);
        for (const size_t neigh : r_adj[idx])
            q.push(neigh);
    }

    // find all the So to keep -> i.e. are target, or is used outside the
    // So_group
    std::vector<size_t> So;
    // note that we construct So from So_group because some entry in So_group
    // might be used by some addition outside So_group, so we need it to still
    // be an intermediate sum
    for (const size_t &idx : So_group) {
        // we don't want the inputs to be part of the outputs
        if (idx < G.n)
            continue;

        // keep target nodes
        if (output_idxs.count(idx)) {
            So.push_back(idx);
            continue;
        }

        // keep nodes that are used outside So_group
        std::vector<size_t> &built_on_idx_vec = adj[idx];
        for (size_t &built_on_idx : built_on_idx_vec) {
            if (So_group.count(built_on_idx))
                continue;
            So.push_back(idx);
            break;
        }
    }
    std::sort(So.begin(), So.end());

    // construct the input
    std::vector<size_t> Si;
    assert(q.empty()); // sanity check
    for (size_t idx : So)
        q.push(idx);

    std::unordered_set<size_t> seen;
    while (!q.empty()) {
        size_t idx = q.front();
        q.pop();
        if (seen.count(idx))
            continue;
        seen.insert(idx);
        if (idx < G.n)
            Si.push_back(idx);
        for (size_t neigh : r_adj[idx])
            q.push(neigh);
    }
    std::sort(Si.begin(), Si.end());

    // construct the new G, rows are targets (So) made from Si
    size_t m = So.size();
    size_t n = Si.size();
    assert(m <= 64 && n <= 64); // sanity check
    std::vector<uint64_t> matrix(n, 0);
    for (size_t col = 0; col < n; col++)
        for (size_t row = 0; row < m; row++)
            if (basis[So[row]] & basis[Si[col]])
                matrix[col] |= 1ULL << row; // each entry of matrix is a column

    return {Z2Matrix(matrix, m, n), Si, So};
}

Result merge_results(const Z2Matrix &G, const Result &result_G,
                     const Z2Matrix &new_G, const Result &result_new_G,
                     const std::vector<size_t> &Si,
                     const std::vector<size_t> &So) {
    assert(result_new_G.method.outputs.size() == So.size());

    Result result;
    result.additions_before = result_G.additions_before;

    // setup So_group
    std::vector<std::vector<size_t>> r_adj(G.n +
                                           result_G.method.additions.size());
    for (size_t i = 0; i < result_G.method.additions.size(); i++) {
        size_t idx = G.n + i;
        const auto &[idx1, idx2] = result_G.method.additions[i];
        r_adj[idx].push_back(idx1);
        r_adj[idx].push_back(idx2);
    }
    std::unordered_set<size_t> So_group;
    std::queue<size_t> q;
    for (const size_t &idx : So)
        q.push(idx);
    while (!q.empty()) {
        size_t idx = q.front();
        q.pop();
        if (So_group.count(idx))
            continue;
        So_group.insert(idx);
        for (const size_t neigh : r_adj[idx])
            q.push(neigh);
    }

    // convert result_new_G back into n-dimensional space
    for (size_t i = 0; i < result_new_G.method.additions.size(); i++) {
        auto &[_idx1, _idx2] = result_new_G.method.additions[i];
        size_t idx1 = _idx1 < new_G.n ? Si[_idx1] : G.n - new_G.n + _idx1;
        size_t idx2 = _idx2 < new_G.n ? Si[_idx2] : G.n - new_G.n + _idx2;
        result.method.additions.push_back({idx1, idx2});
    }

    // small helper for mapping G indices to updated new_G indices i.e.
    // after new_G has been expanded back into G.n dimensions
    std::unordered_map<size_t, size_t> old_idx2new_idx;
    for (size_t i = 0; i < result_new_G.method.outputs.size(); i++) {
        size_t o = result_new_G.method.outputs[i];
        old_idx2new_idx[So[i]] = o < new_G.n ? Si[o] : G.n - new_G.n + o;
    }

    // merge the two addition methods
    for (size_t i = 0; i < result_G.method.additions.size(); i++) {
        size_t idx = i + G.n;
        if (So_group.count(idx))
            continue;
        const auto &[_idx1, _idx2] = result_G.method.additions[i];
        size_t idx1 = _idx1, idx2 = _idx2;
        std::array<std::reference_wrapper<size_t>, 2> c_idxs = {idx1, idx2};
        for (size_t &c_idx : c_idxs) {
            if (So_group.count(c_idx)) {
                // sanity check
                assert(old_idx2new_idx.count(c_idx));
                // c_idx in So, we map to new index (output in new_G)
                c_idx = old_idx2new_idx[c_idx];
            } else if (c_idx < G.n) {
                // don't do anything to indices in E
                continue;
            } else {
                // remap idxs according to computing So part first
                size_t num_old_So_group_additions = So_group.size() - Si.size();
                size_t num_new_So_group_additions =
                    result_new_G.method.additions.size();
                c_idx = (c_idx + num_new_So_group_additions) -
                        num_old_So_group_additions;
            }
        }
        result.method.additions.push_back({idx1, idx2});
    }
    result.additions_after = result.method.additions.size();

    // fill the output vector
    for (const size_t &o : result_G.method.outputs) {
        if (o < G.n || o == std::numeric_limits<size_t>::max()) {
            result.method.outputs.push_back(o);
        } else if (So_group.count(o)) {
            result.method.outputs.push_back(old_idx2new_idx[o]);
        } else {
            result.method.outputs.push_back(o - (G.n - new_G.n));
        }
    }

    // finally, topo-sort the additions
    toposorter(result.method, G.n);

    return result;
}
} // namespace slp::gf2
