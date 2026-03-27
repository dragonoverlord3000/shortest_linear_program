#include "slp/boyar_peralta/internal.hpp"
#include "slp/types.hpp"

#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// this implementation is inspired by:
// https://github.com/rub-hgi/shorter_linear_slps_for_mds_matrices/blob/master/slp_heuristic.cpp
// https://eprint.iacr.org/2025/1493.pdf

namespace slp::gf2 {

namespace {

class Basis {
  private:
    // TODO: find a data-based threshold, possibly based also on density and
    // such
    slp::ReachableStrategy reachable_strategy;

  public:
    std::size_t n; // the dimension size
    std::unordered_set<uint64_t> s_basis;
    std::vector<uint64_t> basis;

    Basis(std::size_t n, ReachableStrategy reachable_strategy)
        : reachable_strategy(reachable_strategy), n(n) {}

    void add_element(uint64_t b) {
        s_basis.insert(b);
        basis.push_back(b);
    }

    std::size_t get_dist(uint64_t t, uint64_t new_b, std::size_t prev_dist) {
        assert(!basis.empty()); // otherwise basis.size() - 1 gives something
                                // completely wrong

        uint64_t new_target = t ^ new_b;
        switch (prev_dist) {
        case 0:
            return 0;
        case 1:
            return new_target == 0 ? 0 : 1;
        case 2:
            return s_basis.count(new_target) ? 1 : 2;
        case 3:
            for (uint64_t b : basis)
                if (s_basis.count(new_target ^ b))
                    return 2;
            return 3;
        }

        // TODO: implement a fast (ISD) method like Stern/Dumer
        switch (reachable_strategy) {
        case slp::ReachableStrategy::BruteForce:
            return _brute_reachable(t ^ new_b, basis.size() - 1, prev_dist - 1)
                       ? prev_dist - 1
                       : prev_dist;
        case slp::ReachableStrategy::MITM:
            return _mitm_reachable(t ^ new_b, prev_dist - 1) ? prev_dist - 1
                                                             : prev_dist;
        case slp::ReachableStrategy::BacktracingSparseAware: {
            std::vector<std::unordered_set<std::size_t>> column2setbasisidxs(n);
            for (std::size_t basis_idx = 0; basis_idx < basis.size();
                 basis_idx++)
                for (std::size_t shift = 0; shift < n; shift++)
                    if (basis[basis_idx] & (1ULL << shift))
                        column2setbasisidxs[shift].insert(basis_idx);

            return _sparse_aware_bt_reachable(t ^ new_b, prev_dist - 1,
                                              column2setbasisidxs)
                       ? prev_dist - 1
                       : prev_dist;
        }
        default:
            assert(false); // not a valid case
        }

        return prev_dist;
    }

    bool _sparse_aware_bt_reachable(
        uint64_t t, std::size_t dist,
        std::vector<std::unordered_set<std::size_t>> &column2setbasisidxs) {
        if (t == 0)
            return true;
        if (dist == 0)
            return false;
        if (dist == 1)
            return s_basis.count(
                t); // note that it is okay to reuse basis elements, but it will
                    // never be what is required (by construction)

        std::size_t best_idx = n;
        for (std::size_t i = 0; i < n; i++) {
            if ((t & (1ULL << i)) == 0)
                continue;
            if (!column2setbasisidxs[i].empty() &&
                (best_idx == n || column2setbasisidxs[i].size() <
                                      column2setbasisidxs[best_idx].size()))
                best_idx = i;
        }

        if (best_idx == n)
            return false;

        std::vector<std::size_t> basis_idxs(
            column2setbasisidxs[best_idx].begin(),
            column2setbasisidxs[best_idx].end());

        assert(basis_idxs.size() <
               64); // if this assert ever fails, we shouldn't be using
                    // _sparse_aware_bt_reachable anyway for this case
        for (uint64_t bm = 1; bm < (1ULL << basis_idxs.size()); bm++) {
            std::size_t pc = std::popcount(bm);
            if (pc % 2 == 0 || pc > dist)
                continue; // TODO: can probably speed this up later by computing
                          // the odd-sized combinations directly
            std::vector<std::size_t> subset;
            uint64_t b = 0;
            for (uint64_t shift = 0; shift < basis_idxs.size(); shift++) {
                if (!(bm & (1ULL << shift)))
                    continue;
                std::size_t idx = basis_idxs[shift];
                subset.push_back(idx);
                b ^= basis[idx];
            }

            // remove used basis vectors, save where removed so as to backtrack
            // later
            std::vector<std::vector<std::size_t>> memory(n);
            for (std::size_t i = 0; i < n; i++)
                for (std::size_t idx : subset)
                    if (column2setbasisidxs[i].count(idx)) {
                        memory[i].push_back(idx);
                        column2setbasisidxs[i].erase(idx);
                    }

            // recurse
            bool ok = _sparse_aware_bt_reachable(t ^ b, dist - pc,
                                                 column2setbasisidxs);
            if (ok)
                return true;

            // backtrack
            for (std::size_t i = 0; i < n; i++)
                for (std::size_t idx : memory[i])
                    column2setbasisidxs[i].insert(idx);
        }

        return false;
    }

    bool _mitm_reachable(uint64_t t, std::size_t dist) {
        std::unordered_map<uint64_t, std::size_t> even{{0, 0}}, odd{{0, 0}};
        for (std::size_t i = 0; i < basis.size(); i++) {
            std::unordered_map<uint64_t, std::size_t> &s = (i & 1) ? odd : even;
            std::vector<std::pair<uint64_t, std::size_t>> updates;
            for (auto &[b, d] : s) {
                if (d >= dist)
                    continue;
                updates.push_back({b ^ basis[i], d + 1});
                if ((b ^ basis[i]) == t)
                    return true;
            }
            for (auto &[b, d] : updates) {
                if (s.count(b) && s[b] <= d)
                    continue;
                s[b] = d;
            }
        }

        // note |odd| <= |even| as we index starting from 0
        for (auto &[b, d] : odd) {
            uint64_t new_t = b ^ t;
            if (!even.count(new_t))
                continue;
            if (even[new_t] + d <= dist)
                return true;
        }

        return false;
    }

    bool _brute_reachable(uint64_t t, std::size_t at, std::size_t d) const {
        if (d == 0)
            return t == 0;
        if (at + 1 < d)
            return false;
        if (d == 1) {
            for (std::size_t i = 0; i <= at; i++)
                if (basis[i] == t)
                    return true;
            return false;
        }

        if (_brute_reachable(t ^ basis[at], at - 1, d - 1) ||
            _brute_reachable(t, at - 1, d))
            return true;
        return false;
    }

    bool contains(uint64_t b) const { return s_basis.count(b); }
    std::size_t size() const { return basis.size(); }

    uint64_t operator[](std::size_t idx) const { return basis[idx]; }
};

void apply_move(Basis &basis, std::vector<std::size_t> &new_dist,
                std::vector<std::size_t> &dist,
                std::vector<std::pair<std::size_t, std::size_t>> &additions,
                std::size_t i, std::size_t j, uint64_t new_b) {
    dist = new_dist;
    additions.push_back({i, j});
    basis.add_element(new_b);
}

std::pair<std::size_t, std::size_t>
evaluate_move(Basis &basis, const std::vector<uint64_t> &targets,
              std::vector<std::size_t> &new_dist,
              const std::vector<std::size_t> &prev_dist, const uint64_t new_b) {
    std::size_t cur_d = 0, cur_nd = 0;
    new_dist.assign(prev_dist.size(), 0);

    for (std::size_t idx = 0; idx < targets.size(); idx++) {
        if (new_b == targets[idx] || prev_dist[idx] == 0)
            continue;

        std::size_t d = basis.get_dist(targets[idx], new_b, prev_dist[idx]);
        new_dist[idx] = d;
        cur_d += new_dist[idx];
        cur_nd += new_dist[idx] * new_dist[idx];
    }

    return {cur_d, cur_nd};
}

void step(Basis &basis, const std::vector<uint64_t> &targets,
          std::vector<std::size_t> &dist, std::size_t m,
          std::vector<std::pair<std::size_t, std::size_t>> &additions,
          std::unordered_set<uint64_t> &s_targets_missing) {
    std::size_t best_dist_norm = 0;
    std::size_t best_dist_sum = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> best_dist(m);
    std::size_t best_i = 0;
    std::size_t best_j = 0;

    for (std::size_t i = 0; i < basis.size(); i++) {
        for (std::size_t j = i + 1; j < basis.size(); j++) {
            uint64_t new_b = basis[i] ^ basis[j];
            if (basis.contains(new_b))
                continue;
            // if dist[some_target] = 1 then make dist[some_target] = 0
            if (s_targets_missing.count(new_b)) {
                std::vector<std::size_t> new_dist;
                evaluate_move(basis, targets, new_dist, dist, new_b);
                apply_move(basis, new_dist, dist, additions, i, j, new_b);
                s_targets_missing.erase(new_b);
                return;
            }

            std::vector<std::size_t> new_dist;
            auto [cur_d, cur_nd] =
                evaluate_move(basis, targets, new_dist, dist, new_b);
            if ((cur_d < best_dist_sum) ||
                (cur_d == best_dist_sum && cur_nd > best_dist_norm)) {
                best_dist_sum = cur_d;
                best_dist_norm = cur_nd;
                best_dist = new_dist;
                best_i = i;
                best_j = j;
            }
        }
    }

    uint64_t best_b = basis[best_i] ^ basis[best_j];
    apply_move(basis, best_dist, dist, additions, best_i, best_j, best_b);
}
} // namespace

// returns method
// pair `p_i` in method means that element `i` is constructed by taking
// B[p_i[0]] xor B[p_i[1]], where B is the basis
std::vector<std::pair<std::size_t, std::size_t>>
run_boyar_peralta(const std::vector<uint64_t> &G, std::size_t m, std::size_t n,
                  const slp::Options &options) {
    assert(m <= 64 && n <= 64);

    // each row of G is a target, G[j] is a column (i.e. variable)
    Basis basis(n, options.reachable_strategy);
    std::unordered_set<uint64_t> s_targets_missing;
    for (std::size_t shift = 0; shift < n; shift++) {
        uint64_t b = 1ULL << shift;
        basis.add_element(b);
    }

    std::vector<uint64_t> targets(m);
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < n; j++) {
            if (G[j] & (1ULL << i))
                targets[i] |= 1ULL << j;
        }
        if (targets[i] == 0 || basis.contains(targets[i]))
            continue;
        s_targets_missing.insert(targets[i]);
    }

    std::vector<std::size_t> dist(m);
    for (std::size_t idx = 0; idx < m; idx++) {
        if (targets[idx] == 0 || basis.contains(targets[idx])) {
            dist[idx] = 0;
            continue;
        }
        std::size_t d = 0;
        for (std::size_t col = 0; col < n; col++)
            d += (G[col] & (1ULL << idx)) > 0;
        dist[idx] = d - 1;
    }

    std::vector<std::pair<std::size_t, std::size_t>> additions;
    int num_rounds = 0;
    while (!s_targets_missing.empty()) {
        num_rounds++;
        if (options.verbose)
            std::cout << "BP Round #" << num_rounds << std::endl;

        step(basis, targets, dist, m, additions, s_targets_missing);
    }

    return additions;
}
} // namespace slp::gf2
