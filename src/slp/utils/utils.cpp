#include "slp/utils/utils.hpp"

#include <queue>
#include <unordered_map>

namespace slp::gf2 {
// changes addition_method inplace so that additions are topological sorting
// outputs are updated accordingly
void toposorter(AdditionMethod &addition_method, const size_t n) {
    size_t N = n + addition_method.additions.size();
    std::vector<std::vector<size_t>> adj(N);
    std::vector<size_t> in_degree(N);
    for (size_t idx = 0; idx < addition_method.additions.size(); idx++) {
        auto [idx1, idx2] = addition_method.additions[idx];
        adj[idx1].push_back(idx + n);
        adj[idx2].push_back(idx + n);
        in_degree[idx + n] += 2;
    }

    std::queue<size_t> q;
    for (size_t node = 0; node < N; node++)
        if (in_degree[node] == 0)
            q.push(node);

    std::vector<size_t> new_ordering;
    while (!q.empty()) {
        size_t idx = q.front();
        q.pop();
        if (idx >= n)
            new_ordering.push_back(idx - n);
        for (size_t neigh : adj[idx]) {
            in_degree[neigh]--;
            if (in_degree[neigh] == 0)
                q.push(neigh);
        }
    }

    // sanity check
    assert(new_ordering.size() == addition_method.additions.size());

    // permute additions and outputs according to toposort
    std::unordered_map<size_t, size_t> old2new;
    std::vector<std::pair<size_t, size_t>> additions(
        addition_method.additions.size());
    for (size_t i = 0; i < addition_method.additions.size(); i++) {
        additions[i] = addition_method.additions[new_ordering[i]];
        old2new[new_ordering[i] + n] = i + n;
    }
    for (size_t i = 0; i < addition_method.additions.size(); i++) {
        auto &[idx1, idx2] = additions[i];
        if (idx1 >= n)
            idx1 = old2new[idx1];
        if (idx2 >= n)
            idx2 = old2new[idx2];
        if (idx2 < idx1)
            std::swap(idx1,
                      idx2); // just nice aesthetically to have them ordered
    }
    // overwrite with the new additions
    addition_method.additions = std::move(additions);

    for (size_t i = 0; i < addition_method.outputs.size(); i++) {
        size_t &o = addition_method.outputs[i];
        // SIZE_T_MAX denotes a zero-row
        if (o == std::numeric_limits<size_t>::max() || o < n)
            continue;
        o = old2new[o];
    }
}
} // namespace slp::gf2
