#include "slp/framework/framework.hpp"
#include "slp/types.hpp"
#include "slp/utils/utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>

namespace slp::gf2::fw {

namespace {

enum class VisitState { Visiting, Done };

void fill_basis2idx_checked(
    uint64_t target,
    const std::unordered_map<uint64_t, std::vector<uint64_t>> &B_adj,
    std::unordered_map<uint64_t, size_t> &basis2idx,
    std::vector<std::pair<size_t, size_t>> &additions, size_t n,
    std::unordered_map<uint64_t, VisitState> &state) {
    if (basis2idx.count(target))
        return;

    if (state.count(target) && state[target] == VisitState::Visiting) {
        std::cerr << "Cycle in B_adj at basis value " << target << std::endl;
        std::abort();
    }

    auto it = B_adj.find(target);
    if (it == B_adj.end()) {
        std::cerr << "Missing B_adj entry for basis value " << target
                  << std::endl;
        std::abort();
    }

    if (it->second.size() != 2) {
        std::cerr << "Bad B_adj arity for basis value " << target
                  << ": expected 2, got " << it->second.size() << std::endl;
        std::abort();
    }

    state[target] = VisitState::Visiting;

    std::vector<size_t> atom_idxs;
    uint64_t checker = 0;

    for (const auto &from_b : it->second) {
        fill_basis2idx_checked(from_b, B_adj, basis2idx, additions, n, state);
        checker ^= from_b;
        atom_idxs.push_back(basis2idx.at(from_b));
    }

    if (checker != target) {
        std::cerr << "Bad B_adj equation for target " << target << std::endl;
        std::abort();
    }

    size_t idx = basis2idx.size();

    if (atom_idxs[0] > atom_idxs[1])
        std::swap(atom_idxs[0], atom_idxs[1]);

    additions.push_back({atom_idxs[0], atom_idxs[1]});
    basis2idx[target] = idx;

    state[target] = VisitState::Done;
}
} // namespace

std::tuple<Z2Matrix, std::vector<size_t>, std::vector<size_t>>
construct_new_G(const Z2Matrix &G, const Result &result, const size_t gap,
                const size_t start, const Options &options) {
    if (options.debug)
        validate_method(result.method, G.n, "construct_new_G input");
    assert(G.n <= 64 && G.m <= 64);

    // construct the unprocessed set So
    std::unordered_set<size_t> So_set;
    for (size_t idx = start; idx < start + gap; idx++)
        So_set.insert(idx + G.n);

    // remove elements that are not used outside So and are not targets
    std::unordered_set<size_t> target_set(result.method.outputs.begin(),
                                          result.method.outputs.end()),
        used_outside_set;
    for (size_t i = 0; i < result.method.additions.size(); i++) {
        size_t idx = i + G.n;
        auto &[idx1, idx2] = result.method.additions[i];
        if (!So_set.count(idx)) {
            used_outside_set.insert(idx1);
            used_outside_set.insert(idx2);
        }
    }

    std::vector<size_t> So;
    for (const size_t idx : So_set) {
        if (target_set.count(idx) || used_outside_set.count(idx))
            So.push_back(idx);
    }
    std::sort(So.begin(), So.end());
    if (So.empty()) {
        std::cerr << "construct_new_G produced empty So" << ", gap=" << gap
                  << ", start=" << start << std::endl;
    }

    // find Si
    std::unordered_map<size_t, std::unordered_set<size_t>> So_sums(So.size());
    std::unordered_set<size_t> Si_set;
    for (size_t idx : So) {
        std::vector<size_t> stack = {idx};
        while (!stack.empty()) {
            size_t node = stack.back();

            if (options.debug) {
                if (node >= G.n + result.method.additions.size()) {
                    std::cerr
                        << "Invalid node in construct_new_G: node=" << node
                        << ", total_basis_size="
                        << (G.n + result.method.additions.size())
                        << ", gap=" << gap << ", start=" << start << std::endl;
                    std::abort();
                }
            }

            stack.pop_back();
            if (So_set.count(node)) {
                stack.push_back(result.method.additions[node - G.n].first);
                stack.push_back(result.method.additions[node - G.n].second);
            } else {
                Si_set.insert(node);
                if (So_sums[idx].count(node))
                    So_sums[idx].erase(node); // a ^ a = 0
                else
                    So_sums[idx].insert(node);
            }
        }
    }
    /*
    // how the original framework paper does it, seems to work worse
    for(size_t j = 0; j < So.size(); j++)
        for(size_t i = 0; i < G.n; i++)
            if (So[j] & (1ULL << i)) {
                std::unordered_set<size_t>& So_sum = So_sums[So[j]];
                So_sum.insert(i);
                Si_set.insert(i);
            }
            */

    std::vector<size_t> Si(Si_set.begin(), Si_set.end());
    std::sort(Si.begin(), Si.end());

    // Construct G
    size_t m = So.size();
    size_t n = Si.size();
    if (m > 64 || n > 64) {
        std::vector<size_t> basis_Si;
        for (size_t i = 0; i < G.n; i++)
            basis_Si.push_back(i);
        std::vector<size_t> basis_So;
        for (const size_t o : result.method.outputs)
            basis_So.push_back(o);
        return {G, basis_Si, basis_So};
    }
    assert(m <= 64 && n <= 64);

    std::vector<uint64_t> G_prime(n, 0);
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++) {
            if (So_sums[So[i]].count(Si[j]))
                G_prime[j] |= 1ULL << i;
        }

    return {Z2Matrix(G_prime, m, n), Si, So};
}

Result merge_results(const Z2Matrix &G, const Result &result_G,
                     const Z2Matrix &new_G, const Result &result_new_G,
                     const std::vector<size_t> &Si,
                     const std::vector<size_t> &So, const Options &options) {
    assert(result_new_G.method.outputs.size() == So.size());
    assert(G.m <= 64 && G.n <= 64);
    assert(new_G.m <= 64 && new_G.n <= 64);
    assert(Si.size() == new_G.n);
    assert(So.size() == new_G.m);
    validate_method(result_G.method, G.n, "merge_results result_G input");
    validate_method(result_new_G.method, new_G.n,
                    "merge_results result_new_G input");

    std::vector<std::vector<size_t>> adj(G.n +
                                         result_G.method.additions.size());
    for (size_t i = 0; i < result_G.method.additions.size(); i++) {
        size_t idx = i + G.n;
        auto &[idx1, idx2] = result_G.method.additions[i];
        adj.at(idx1).push_back(idx);
        adj.at(idx2).push_back(idx);
    }
    std::unordered_set<size_t> So_set(So.begin(), So.end()),
        target_set(result_G.method.outputs.begin(),
                   result_G.method.outputs.end());

    // form the basis for G and new_G
    std::vector<uint64_t> basis_G(G.n), basis_new_G(new_G.n);
    for (size_t i = 0; i < G.n; i++)
        basis_G[i] |= 1ULL << i;
    for (auto &[idx1, idx2] : result_G.method.additions)
        basis_G.push_back(basis_G[idx1] ^ basis_G[idx2]);
    for (size_t i = 0; i < new_G.n; i++)
        basis_new_G[i] = basis_G[Si[i]];
    for (auto &[idx1, idx2] : result_new_G.method.additions)
        basis_new_G.push_back(basis_new_G[idx1] ^ basis_new_G[idx2]);

    // form the dependency graph based on the values of B, rather than the
    // indices, I.e. basis[idx] = basis[idx1] ^ basis[idx2]
    std::unordered_map<uint64_t, std::vector<uint64_t>> B_adj;

    // first we do the part of G that builds Si
    std::vector<size_t> Si_stack(Si.begin(), Si.end());
    while (!Si_stack.empty()) {
        size_t idx = Si_stack.back();
        Si_stack.pop_back();

        uint64_t b = basis_G[idx];
        // skip values we already have calculated as part of computing Si
        if (B_adj.count(b))
            continue;

        // idx is one of the original input nodes
        if (idx < G.n)
            continue;
        auto &[idx1, idx2] = result_G.method.additions[idx - G.n];

        // sanity check
        assert(b == (basis_G[idx1] ^ basis_G[idx2]));

        Si_stack.push_back(idx1);
        Si_stack.push_back(idx2);

        B_adj[b].push_back(basis_G[idx1]);
        B_adj[b].push_back(basis_G[idx2]);
    }

    // now we do new_G, because we know we want to keep these
    for (size_t i = 0; i < result_new_G.method.additions.size(); i++) {
        size_t idx = i + new_G.n;
        uint64_t b = basis_new_G[idx];

        // skip values we already have as part of new_G or part of computing Si
        if (B_adj.count(b))
            continue;

        auto &[idx1, idx2] = result_new_G.method.additions[i];

        // sanity check
        assert(b == (basis_new_G[idx1] ^ basis_new_G[idx2]));
        // std::cout << b << "| " << basis_new_G[idx1] << "| " <<
        // basis_new_G[idx2]
        //           << std::endl;

        B_adj[b].push_back(basis_new_G[idx1]);
        B_adj[b].push_back(basis_new_G[idx2]);
    }

    for (size_t i = 0; i < result_G.method.additions.size(); i++) {
        size_t idx = i + G.n;
        uint64_t b = basis_G[idx];

        // this step also ensure that we avoid cycles
        if (B_adj.count(b))
            continue; // already covered (by new_G)

        auto &[idx1, idx2] = result_G.method.additions[i];

        // sanity check
        assert(b == (basis_G[idx1] ^ basis_G[idx2]));
        // std::cout << b << ", " << basis_G[idx1] << ", " << basis_G[idx2]
        //           << std::endl;

        B_adj[b].push_back(basis_G[idx1]);
        B_adj[b].push_back(basis_G[idx2]);
    }
    // whenever we describe something using new_G, remove all old_G descriptions
    // this is probably redundant
    if (options.debug) {
        for (auto &[_, v] : B_adj) {
            // sanity check
            assert(v.size() == 2 && "B_adj creates key with not 2 values");
        }
    }

    // finally, reorder back into a result format
    Result result;
    result.additions_before = result_G.additions_before;

    std::unordered_map<uint64_t, size_t> basis2idx;
    for (size_t i = 0; i < G.n; i++)
        basis2idx[1ULL << i] = i;

    std::unordered_map<uint64_t, VisitState> state;
    for (const size_t target_idx : target_set) {
        if (target_idx == std::numeric_limits<size_t>::max())
            continue; // skip zero row targets
        uint64_t target = basis_G[target_idx];
        fill_basis2idx_checked(target, B_adj, basis2idx,
                               result.method.additions, G.n, state);
    }

    // note that it doesn't matter what the value of result_new_G.method.outputs
    // is, since we're going solely of off the basis value graph in stead
    for (const size_t o : result_G.method.outputs) {
        if (o == std::numeric_limits<size_t>::max())
            result.method.outputs.push_back(o);
        else
            result.method.outputs.push_back(basis2idx[basis_G[o]]);
    }

    result.additions_after = result.method.additions.size();

    // its nice to have them in topological order
    toposorter(result.method, G.n);

    return result;
}
} // namespace slp::gf2::fw
