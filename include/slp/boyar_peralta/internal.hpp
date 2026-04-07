#pragma once

#include "slp/types.hpp"

#include <unordered_map>
#include <unordered_set>

namespace slp::gf2::bp {
class Basis {
  private:
    // TODO: find a data-based threshold, possibly based also on density and
    // such
    slp::ReachableStrategy reachable_strategy;

  public:
    size_t n; // the dimension size
    std::unordered_set<uint64_t> s_basis;
    std::vector<uint64_t> basis;
    std::vector<std::unordered_set<size_t>> column2setbasisidxs;

    Basis(size_t n, ReachableStrategy reachable_strategy)
        : reachable_strategy(reachable_strategy), n(n) {
        if (reachable_strategy ==
            slp::ReachableStrategy::BacktracingSparseAware)
            column2setbasisidxs.assign(n, {});
    }

    void add_element(uint64_t b) {
        s_basis.insert(b);
        basis.push_back(b);

        // id we're using sparse aware backtracking then update the
        // column2setbasisidxs when adding elements
        if (reachable_strategy ==
            slp::ReachableStrategy::BacktracingSparseAware) {
            uint64_t x = b;
            size_t idx = basis.size() - 1;
            while (x) {
                size_t tz = std::countr_zero(x);
                column2setbasisidxs[tz].insert(idx);
                x &= x - 1;
            }
        }
    }

    size_t get_dist(uint64_t t, uint64_t new_b, size_t prev_dist) {
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
        case slp::ReachableStrategy::BacktracingSparseAware:
            return _sparse_aware_bt_reachable(t ^ new_b, prev_dist - 1,
                                              column2setbasisidxs)
                       ? prev_dist - 1
                       : prev_dist;
        default:
            assert(false); // not a valid case
        }

        return prev_dist;
    }

    bool _sparse_aware_bt_reachable(
        uint64_t t, size_t dist,
        std::vector<std::unordered_set<size_t>> &column2setbasisidxs) {
        if (t == 0)
            return true;
        if (dist == 0)
            return false;
        if (dist == 1)
            return s_basis.count(
                t); // note that it is okay to reuse basis elements, but it will
                    // never be what is required (by construction)

        size_t best_idx = n;
        uint64_t x = t;
        while (x) {
            size_t tz = std::countr_zero(x);

            // t has a bit at idx i, but no way to set it with any basis
            if (column2setbasisidxs[tz].empty())
                return false;

            // if no best_idx found so far, or i is better idx than current
            // `best_idx`
            if (best_idx == n || column2setbasisidxs[tz].size() <
                                     column2setbasisidxs[best_idx].size())
                best_idx = tz;

            x &= x - 1; // remove the least significant set bit trick
        }

        std::vector<size_t> basis_idxs(column2setbasisidxs[best_idx].begin(),
                                       column2setbasisidxs[best_idx].end());
        // heuristically, we want larger overlaps with t to be first -> didn't
        // seem to have any tangible effect when benchmarking
        // std::sort(basis_idxs.begin(), basis_idxs.end(), [&t](const size_t& a,
        // const size_t& b) { return std::popcount(a ^ t) < std::popcount(b ^
        // t); });

        assert(basis_idxs.size() <
               64); // if this assert ever fails, we shouldn't be using
                    // _sparse_aware_bt_reachable anyway for this case
        for (uint64_t bm = 1; bm < (1ULL << basis_idxs.size()); bm++) {
            size_t pc = std::popcount(bm);
            if (pc % 2 == 0 || pc > dist)
                continue; // TODO: can probably speed this up later by computing
                          // the odd-sized combinations directly, and starting
                          // from smallest number of set bits 1, then moving up
                          // to 3 then after trying all those moving up to 5,
                          // ...
            std::vector<size_t> subset;
            uint64_t b = 0;
            uint64_t x = bm;
            while (x) {
                size_t idx = basis_idxs[std::countr_zero(x)];
                subset.push_back(idx);
                b ^= basis[idx];
                x &= x - 1;
            }

            // remove used basis vectors, save where removed so as to backtrack
            // later
            std::vector<std::vector<size_t>> memory(n);
            for (size_t i = 0; i < n; i++)
                for (size_t idx : subset)
                    if (column2setbasisidxs[i].count(idx)) {
                        memory[i].push_back(idx);
                        column2setbasisidxs[i].erase(idx);
                    }

            // recurse
            bool ok = _sparse_aware_bt_reachable(t ^ b, dist - pc,
                                                 column2setbasisidxs);
            // backtrack
            for (size_t i = 0; i < n; i++)
                for (size_t idx : memory[i])
                    column2setbasisidxs[i].insert(idx);

            // if worked, then return so, otherwise keep trying other
            // combinations
            if (ok)
                return true;
        }

        return false;
    }

    bool _mitm_reachable(uint64_t t, size_t dist) {
        std::unordered_map<uint64_t, size_t> even{{0, 0}}, odd{{0, 0}};
        for (size_t i = 0; i < basis.size(); i++) {
            std::unordered_map<uint64_t, size_t> &s = (i & 1) ? odd : even;
            std::vector<std::pair<uint64_t, size_t>> updates;
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

    bool _brute_reachable(uint64_t t, size_t at, size_t d) const {
        if (d == 0)
            return t == 0;
        if (at + 1 < d)
            return false;
        if (d == 1) {
            for (size_t i = 0; i <= at; i++)
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
    size_t size() const { return basis.size(); }

    uint64_t operator[](size_t idx) const { return basis[idx]; }
};

void apply_move_bp(Basis &basis, std::vector<size_t> &new_dist,
                   std::vector<size_t> &dist,
                   std::vector<std::pair<size_t, size_t>> &additions, size_t i,
                   size_t j, uint64_t new_b);

std::pair<size_t, size_t> evaluate_move_bp(Basis &basis,
                                           const std::vector<uint64_t> &targets,
                                           std::vector<size_t> &new_dist,
                                           const std::vector<size_t> &prev_dist,
                                           const uint64_t new_b);

bool evaluate_move_Ax_filter(Basis &basis, const std::vector<uint64_t> &targets,
                             std::vector<size_t> &new_dist,
                             const std::vector<size_t> &prev_dist,
                             const uint64_t new_b,
                             const std::vector<size_t> &filter_indices,
                             const bool complement_idxs);

std::pair<size_t, size_t> get_dist_metrics(std::vector<size_t>& dist);

std::vector<std::pair<size_t, size_t>> run_RNBP(const std::vector<uint64_t> &G,
                                                size_t m, size_t n,
                                                const slp::Options &options);

std::vector<std::pair<size_t, size_t>> run_BP(const std::vector<uint64_t> &G,
                                              size_t m, size_t n,
                                              const slp::Options &options);

std::vector<std::pair<size_t, size_t>> run_Ax(const std::vector<uint64_t> &G,
                                              size_t m, size_t n, size_t x,
                                              const slp::Options &options);

AdditionMethod
convert_bp_method(std::vector<uint64_t> &G, size_t m, size_t n,
                  std::vector<std::pair<size_t, size_t>> &additions);

void init_bp(const std::vector<uint64_t> &G, Basis &basis,
             std::unordered_set<uint64_t> &s_targets_missing,
             std::vector<uint64_t> &targets, std::vector<size_t> &dist,
             size_t m, size_t n);
} // namespace slp::gf2::bp
