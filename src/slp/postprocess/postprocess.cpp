#include "slp/postprocess/postprocess.hpp"
#include "slp/utils/utils.hpp"

#include <limits>
#include <unordered_map>
#include <unordered_set>

namespace slp::gf2 {

namespace {

/*
void fill_basis2idx(
const uint64_t target,
const std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>> &B_adj,
std::unordered_map<uint64_t, size_t> &basis2idx,
std::vector<std::pair<size_t, size_t>> &additions, const size_t n) {
if (basis2idx.count(target))
    return;

std::vector<size_t> atom_idxs;
uint64_t checker = 0;

// std::cout << "target: " << target << std::endl;
std::vector<uint64_t> from_bs = {B_adj.at(target).first,
                                 B_adj.at(target).second};
for (auto &from_b : from_bs) {
    fill_basis2idx(from_b, B_adj, basis2idx, additions, n);
    checker ^= from_b;
    atom_idxs.push_back(basis2idx[from_b]);
}

size_t idx = basis2idx.size();
// sanity checks
assert(atom_idxs.size() == 2);
assert(checker == target);
assert(idx == additions.size() + n);

if (atom_idxs[0] > atom_idxs[1])
    std::swap(atom_idxs[0], atom_idxs[1]); // just a nice to have

if (basis2idx.count(target))
    return; // target is already created

// fill in, in post-order traversal manner
additions.push_back({atom_idxs[0], atom_idxs[1]});
basis2idx[target] = idx;
}
*/

// sets up the DAG, removes the temporary variable in basis[idx], tries to
// restore graph
bool temp_var_remover(const std::vector<uint64_t> &basis,
                      const size_t idx_to_remove,
                      const AdditionMethod &addition_method, const size_t n,
                      AdditionMethod &new_addition_method) {
    assert(new_addition_method.additions.empty() &&
           new_addition_method.outputs.empty());
    assert(n <= 64);

    // the indices of the basis vectors to restore
    std::vector<size_t> to_restore;
    // the indices that can be used to restore
    std::vector<size_t> valid_idxs;
    for (size_t i = 0; i < n; i++)
        valid_idxs.push_back(i);
    // the indices we cannot use for restoring, based on the topo ordering
    std::unordered_set<size_t> banned_idxs = {idx_to_remove};
    for (size_t i = 0; i < addition_method.additions.size(); i++) {
        size_t idx = i + n;
        auto &[idx1, idx2] = addition_method.additions[i];
        if (banned_idxs.count(idx) || banned_idxs.count(idx1) ||
            banned_idxs.count(idx2)) {
            if (idx1 == idx_to_remove || idx2 == idx_to_remove)
                to_restore.push_back(idx);
            banned_idxs.insert(idx);
            continue;
        }
        valid_idxs.push_back(idx);
    }

    // see if possible to restore
    std::vector<std::pair<size_t, size_t>> restore_additions;
    for (size_t restore_idx : to_restore) {
        uint64_t s = basis[restore_idx];
        bool found = false;
        size_t restore_i = 0;
        size_t restore_j = 0;
        for (size_t i = 0; i < valid_idxs.size(); i++) {
            for (size_t j = i + 1; j < valid_idxs.size(); j++) {
                if (s != (basis[valid_idxs[i]] ^ basis[valid_idxs[j]]))
                    continue;
                found = true;
                restore_i = valid_idxs[i];
                restore_j = valid_idxs[j];
                break;
            }
            if (found)
                break;
        }
        if (!found)
            return false;
        restore_additions.push_back({restore_i, restore_j});
    }
    // sanity check
    assert(restore_additions.size() == to_restore.size());

    // construct new addition method additions
    size_t restore_pointer = 0;
    for (size_t i = 0; i < addition_method.additions.size(); i++) {
        size_t idx = i + n;
        if (idx == idx_to_remove)
            continue;
        if (restore_pointer < to_restore.size() &&
            to_restore[restore_pointer] == idx) {
            auto [ri, rj] = restore_additions[restore_pointer];
            if (ri > idx_to_remove)
                ri--;
            if (rj > idx_to_remove)
                rj--;
            new_addition_method.additions.push_back({ri, rj});
            restore_pointer++;
            continue;
        }
        auto [idx1, idx2] = addition_method.additions[i];
        if (idx1 > idx_to_remove)
            idx1--;
        if (idx2 > idx_to_remove)
            idx2--;
        new_addition_method.additions.push_back({idx1, idx2});
    }

    // update outputs accordingly
    for (size_t i = 0; i < addition_method.outputs.size(); i++) {
        size_t o = addition_method.outputs[i];
        if (o > idx_to_remove && o != std::numeric_limits<size_t>::max())
            o--;
        new_addition_method.outputs.push_back(o);
    }

    return true;
}

/*
bool temp_var_remover(const std::vector<uint64_t> &_basis,
                      const std::unordered_set<size_t> idxs_to_remove,
                      const AdditionMethod &addition_method, const size_t n,
                      AdditionMethod &new_addition_method) {
    // sanity checks
    assert(new_addition_method.additions.empty() &&
           new_addition_method.outputs.empty());
    assert(idxs_to_remove.size() == 2);
    assert(n <= 64);

    std::vector<uint64_t> basis(_basis.begin(), _basis.end());

    // the indices of the basis vectors to restore
    std::vector<size_t> to_restore;
    // the indices that can be used to restore
    std::vector<size_t> valid_idxs;
    for (size_t i = 0; i < n; i++)
        valid_idxs.push_back(i);
    // the indices we cannot use for restoring, based on the topo ordering
    std::unordered_set<size_t> banned_idxs(idxs_to_remove.begin(),
                                           idxs_to_remove.end());
    for (size_t i = 0; i < addition_method.additions.size(); i++) {
        size_t idx = i + n;
        auto &[idx1, idx2] = addition_method.additions[i];
        if (banned_idxs.count(idx) || banned_idxs.count(idx1) ||
            banned_idxs.count(idx2)) {
            if (idxs_to_remove.count(idx1) || idxs_to_remove.count(idx2))
                to_restore.push_back(idx);
            banned_idxs.insert(idx);
            continue;
        }
        valid_idxs.push_back(idx);
    }

    // see if possible to restore with 1 new addition valid
    bool works = false;
    std::vector<std::pair<size_t, size_t>> restore_additions;
    for (size_t i1 = 0; i1 < valid_idxs.size(); i1++) {
        for (size_t j1 = i1 + 1; j1 < valid_idxs.size(); j1++) {
            uint64_t new_b = basis[valid_idxs[i1]] ^ basis[valid_idxs[j1]];
            size_t extra_idx = basis.size();
            restore_additions.push_back({valid_idxs[i1], valid_idxs[j1]});
            basis.push_back(new_b);
            valid_idxs.push_back(extra_idx);

            works = true;
            for (size_t restore_idx : to_restore) {
                uint64_t s = basis[restore_idx];
                bool found = false;
                size_t restore_i = 0;
                size_t restore_j = 0;
                for (size_t i = 0; i < valid_idxs.size(); i++) {
                    for (size_t j = i + 1; j < valid_idxs.size(); j++) {
                        if (s != (basis[valid_idxs[i]] ^ basis[valid_idxs[j]]))
                            continue;
                        found = true;
                        restore_i = valid_idxs[i];
                        restore_j = valid_idxs[j];
                        break;
                    }
                    if (found)
                        break;
                }
                if (!found) {
                    works = false;
                    break;
                } else
                    restore_additions.push_back({restore_i, restore_j});
            }
            if (works)
                break;

            restore_additions.clear();
            basis.pop_back();
            valid_idxs.pop_back();
        }
        if (works)
            break;
    }
    if (!works)
        return false;

    // construct value-based graph
    std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>> B_adj;
    for (auto &[idx1, idx2] : restore_additions) {
        uint64_t b = basis[idx1] ^ basis[idx2];
        if (B_adj.count(b))
            continue;
        B_adj[b] = {basis[idx1], basis[idx2]};
    }
    for (size_t i = 0; i < addition_method.additions.size(); i++) {
        size_t idx = i + n;
        if (idxs_to_remove.count(idx))
            continue;
        uint64_t b = basis[idx];
        if (B_adj.count(b))
            continue;
        auto &[idx1, idx2] = addition_method.additions[i];
        B_adj[b] = {basis[idx1], basis[idx2]};
    }

    // construct new addition method
    std::unordered_map<uint64_t, size_t> basis2idx;
    for (size_t i = 0; i < n; i++)
        basis2idx[basis[i]] = i;
    for (const size_t o : addition_method.outputs) {
        if (o == std::numeric_limits<size_t>::max()) {
            new_addition_method.outputs.push_back(o);
            continue;
        }
        uint64_t b = basis[o];
        fill_basis2idx(b, B_adj, basis2idx, new_addition_method.additions, n);
        new_addition_method.outputs.push_back(basis2idx[b]);
    }

    return true;
}
*/

} // namespace

AdditionMethod postprocess(const AdditionMethod &addition_method,
                           const size_t n) {
    AdditionMethod processed_addition_method = addition_method;
    bool can_be_removed = true;
    // we start at n bcs removing the inputs doesn't make sense
    while (can_be_removed) {
        can_be_removed = false;

        // ensure basis[i] can be used in forming basis[j] only if i < j,
        // property used in temp_var_remover
        toposorter(processed_addition_method, n);

        // remove single
        for (size_t to_remove = n;
             to_remove < n + processed_addition_method.additions.size();
             to_remove++) {

            // skip removing outputs
            std::unordered_set<size_t> output_set;
            for (const size_t o : processed_addition_method.outputs)
                output_set.insert(o);
            if (output_set.count(to_remove))
                continue;

            std::vector<uint64_t> basis(n);
            for (size_t i = 0; i < n; i++)
                basis[i] = 1ULL << i;
            for (const auto &[idx1, idx2] : processed_addition_method.additions)
                basis.push_back(basis[idx1] ^ basis[idx2]);

            AdditionMethod new_addition_method;
            can_be_removed =
                temp_var_remover(basis, to_remove, processed_addition_method, n,
                                 new_addition_method);
            if (can_be_removed) {
                processed_addition_method = std::move(new_addition_method);
                break;
            }
        }
        if (can_be_removed)
            continue;

        /*
        // ensure basis[i] can be used in forming basis[j] only if i < j,
        // property used in temp_var_remover
        toposorter(processed_addition_method, n);

        // remove double
        for (size_t to_remove_1 = n;
             to_remove_1 < n + processed_addition_method.additions.size();
             to_remove_1++) {
            for (size_t to_remove_2 = to_remove_1 + 1;
                 to_remove_2 < n + processed_addition_method.additions.size();
                 to_remove_2++) {

                // skip removing outputs
                std::unordered_set<size_t> output_set;
                for (const size_t o : processed_addition_method.outputs)
                    output_set.insert(o);
                if (output_set.count(to_remove_1) ||
                    output_set.count(to_remove_2))
                    continue;

                std::vector<uint64_t> basis(n);
                for (size_t i = 0; i < n; i++)
                    basis[i] = 1ULL << i;
                for (const auto &[idx1, idx2] :
                     processed_addition_method.additions)
                    basis.push_back(basis[idx1] ^ basis[idx2]);

                AdditionMethod new_addition_method;
                std::unordered_set<size_t> idxs_to_remove = {to_remove_1,
                                                             to_remove_2};
                can_be_removed = temp_var_remover(basis, idxs_to_remove,
                                                  processed_addition_method, n,
                                                  new_addition_method);
                if (can_be_removed) {
                    processed_addition_method = std::move(new_addition_method);
                    break;
                }
            }
            if (can_be_removed)
                break;
        }
        */
    }
    // not strictly necessary, but it makes the output format nice
    toposorter(processed_addition_method, n);

    return processed_addition_method;
};
} // namespace slp::gf2
