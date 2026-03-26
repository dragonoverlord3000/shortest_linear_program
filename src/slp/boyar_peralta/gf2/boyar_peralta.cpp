#include <slp/boyar_peralta/internal.hpp>

#include <unordered_map>
#include <unordered_set>

// this implementation is inspired by:
// https://github.com/rub-hgi/shorter_linear_slps_for_mds_matrices/blob/master/slp_heuristic.cpp
// https://eprint.iacr.org/2025/1493.pdf

namespace slp::gf2 {

namespace {

class Basis {
  private:
    // every possible construction from the even-index basis vectors,
    // and every possible construction from the odd-index basis vectors.
    // key is # of elements required to construct the values in the value part.
    // Grows in O(2^{n/2}), so for n big enough this will be slow, but still
    // much faster than current methods
    std::unordered_map<std::size_t, std::unordered_set<uint64_t>> even, odd;

  public:
    std::unordered_set<uint64_t> s_basis;
    std::vector<uint64_t> basis;
    Basis() {
        even[0] = {0};
        odd[0] = {0};
    };

    void add_element(uint64_t b) {
        s_basis.insert(b);
        basis.push_back(b);

        std::size_t idx = this->size() - 1;
        std::unordered_map<std::size_t, std::unordered_set<uint64_t>> &s =
            (idx % 2) ? odd : even;
        std::size_t s_size = s.size();
        for (std::size_t sz = s_size; sz > 0; sz--)
            for (uint64_t val : s[sz - 1])
                s[sz].insert(val ^ b);
    }

    // returns min(dist(t) + 1, prev_dist)
    std::size_t get_dist(uint64_t t, std::size_t prev_dist) const {
        // for prev_dist small enough we can do a simple search a bit faster
        if (t == 0)
            return 1;
        if (prev_dist == 2)
            return 2;
        if (prev_dist == 3) {
            if (s_basis.count(t))
                return 2;
            return 3;
        }
        if (prev_dist == 4) {
            for (uint64_t b : basis)
                if (s_basis.count(t ^ b))
                    return 3;
            return 4;
        }

        // otherwise use mitm
        for (std::size_t even_size = 0; even_size < prev_dist - 1;
             even_size++) {
            std::size_t odd_size =
                prev_dist - 2 -
                even_size; // target^new_basis cost 1, even_size + odd_size cost
                           // prev_dist - 2, total cost is prev_dist - 1
            auto ite = even.find(even_size);
            auto ito = odd.find(odd_size);
            if (ite == even.end())
                continue;
            if (ito == odd.end())
                continue;

            for (const uint64_t ov : ito->second)
                if (ite->second.count(t ^ ov))
                    return prev_dist - 1;
        }

        return prev_dist;
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
evaluate_move(const Basis &basis, const std::vector<uint64_t> &targets,
              std::vector<std::size_t> &new_dist,
              const std::vector<std::size_t> &prev_dist, const uint64_t new_b) {
    std::size_t cur_d = 0, cur_nd = 0;
    new_dist.assign(prev_dist.size(), 0);

    for (std::size_t idx = 0; idx < targets.size(); idx++) {
        if (basis.contains(targets[idx]))
            continue;
        std::size_t d = basis.get_dist(new_b ^ targets[idx], prev_dist[idx]);
        new_dist[idx] = d;
        cur_d += new_dist[idx];
        cur_nd += new_dist[idx] * new_dist[idx];
    }

    return {cur_d, cur_nd};
}

void step(Basis &basis, std::vector<uint64_t> &targets,
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
            if (s_targets_missing.count(new_b)) {
                std::vector<std::size_t> new_dist;
                evaluate_move(basis, targets, new_dist, dist, new_b);
                apply_move(basis, new_dist, dist, additions, i, j, new_b);
                s_targets_missing.erase(new_b);
                return;
            }

            std::vector<std::size_t> new_dist(m);
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
run_boyar_peralta(std::vector<uint64_t> &G, std::size_t m, std::size_t n,
                  const slp::Options &options) {

    // each column of G is a target
    Basis basis;
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
    std::vector<std::pair<std::size_t, std::size_t>> method;
    while (!s_targets_missing.empty()) {
        num_rounds++;
        step(basis, targets, dist, m, additions, s_targets_missing);
    }

    return additions;
}
} // namespace slp::gf2
