#include "slp/postprocess/postprocess.hpp"
#include "slp/utils/utils.hpp"

#include <limits>
#include <unordered_set>

namespace slp::gf2 {

namespace {
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
        size_t s = basis[restore_idx];
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

        for (size_t to_remove = n;
             to_remove < n + processed_addition_method.additions.size();
             to_remove++) {
            std::vector<uint64_t> basis(n);
            for (size_t i = 0; i < n; i++)
                basis[i] = 1ULL << i;
            for (const auto &[idx1, idx2] : processed_addition_method.additions)
                basis.push_back(basis[idx1] ^ basis[idx2]);

            std::unordered_set<size_t> output_set;
            for (const size_t o : processed_addition_method.outputs)
                output_set.insert(o);

            if (output_set.count(to_remove))
                continue;
            AdditionMethod new_addition_method;
            can_be_removed =
                temp_var_remover(basis, to_remove, processed_addition_method, n,
                                 new_addition_method);
            if (can_be_removed) {
                processed_addition_method = std::move(new_addition_method);
                break;
            }
        }
    }
    // not strictly necessary, but it makes the output format nice
    toposorter(processed_addition_method, n);

    return processed_addition_method;
};
} // namespace slp::gf2
