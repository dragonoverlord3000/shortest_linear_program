#include "slp/mip/internal.hpp"
#include "slp/potential/internal.hpp"
#include "slp/types.hpp"

#ifdef SLP_WITH_MATHOPT

#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "ortools/math_opt/cpp/math_opt.h"

namespace math_opt = ::operations_research::math_opt;

namespace slp::mip {

Result run_MIP(const Z2Matrix &G, const size_t m, const size_t n) {
    assert(m <= 64 && n < 64);

    std::vector<uint64_t> target_vectors(m);
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            if (G.matrix[j] & (1ULL << i))
                target_vectors[i] |= 1ULL << j;

    math_opt::Model model("mip_model");

    std::map<uint64_t, math_opt::Variable> node_vars;
    const uint64_t N = 1ULL << n;
    for (uint64_t node = 1; node < N; node++)
        node_vars.emplace(
            node, model.AddBinaryVariable(std::format("node_{}", node)));

    std::vector<math_opt::Variable> objective_vec;
    std::map<std::tuple<uint64_t, uint64_t, uint64_t>, math_opt::Variable>
        edges;

    for (uint64_t alpha = 1; alpha < N; alpha++)
        for (uint64_t beta = alpha + 1; beta < N; beta++) {
            const uint64_t gamma = alpha ^ beta;

            if (std::popcount(gamma) <= 1)
                continue;

            const math_opt::Variable t = model.AddBinaryVariable(
                std::format("t_{}_{}_{}", gamma, alpha, beta));

            edges.emplace(
                std::tuple<uint64_t, uint64_t, uint64_t>{gamma, alpha, beta},
                t);
            objective_vec.push_back(t);
        }

    // minimize objective
    math_opt::LinearExpression objective = 0.0;
    for (const math_opt::Variable t : objective_vec)
        objective += t;
    model.Minimize(objective);

    // first set each target vector to 1
    for (uint64_t gamma : target_vectors) {
        if (!node_vars.contains(gamma))
            continue;
        const math_opt::Variable var = node_vars.at(gamma);
        model.AddLinearConstraint(var == 1);
    }

    // set each input vector to 1
    for (size_t i = 0; i < n; i++)
        model.AddLinearConstraint(node_vars.at(1ULL << i) == 1);

    // enforce s_gamma = sum(t_g_a_b)
    for (uint64_t gamma = 1; gamma < N; gamma++) {
        if (std::popcount(gamma) <= 1)
            continue;
        math_opt::LinearExpression incoming = 0.0;
        for (const auto &[key, var] : edges) {
            const uint64_t edge_gamma = std::get<0>(key);
            if (edge_gamma != gamma)
                continue;
            incoming += var;
        }
        model.AddLinearConstraint(node_vars.at(gamma) == incoming);
    }

    // enforce t_g_a_b <= s_a, t_g_a_b <= s_b
    for (const auto &[key, var] : edges) {
        const uint64_t alpha = std::get<1>(key);
        const uint64_t beta = std::get<2>(key);

        model.AddLinearConstraint(var <= node_vars.at(alpha));
        model.AddLinearConstraint(var <= node_vars.at(beta));
    }

    // solve
    math_opt::SolveArguments args;
    args.parameters.enable_output = true;

    const absl::StatusOr<math_opt::SolveResult> solve_result =
        math_opt::Solve(model, math_opt::SolverType::kGscip, args);

    if (!solve_result.ok()) {
        throw std::runtime_error("MathOpt solve failed");
    }

    const absl::Status feasible_status =
        solve_result->termination.EnsureIsOptimalOrFeasible();

    if (!feasible_status.ok()) {
        throw std::runtime_error(
            "MathOpt did not find feasible/optimal solution");
    }

    const auto values = solve_result->variable_values();

    // convert to result format
    Result result;
    std::unordered_set<uint64_t> B;
    for (uint64_t b = 1; b < N; b++) {
        if (values.at(node_vars.at(b)) < 0.5 || B.count(b))
            continue;
        B.insert(b);
    }

    std::vector<uint64_t> full_basis;
    for (size_t i = 0; i < n; i++) {
        full_basis.push_back(1ULL << i);
        B.erase(1ULL << i);
    }

    size_t prev_size = 0;
    while (full_basis.size() > prev_size) {
        prev_size = full_basis.size();
        for (size_t i = 0; i < prev_size; i++)
            for (size_t j = i + 1; j < prev_size; j++) {
                uint64_t b = full_basis[i] ^ full_basis[j];
                if (!B.count(b))
                    continue;
                B.erase(b);
                result.method.additions.push_back({i, j});
                full_basis.push_back(b);
            }
    }

    result.method.outputs.assign(m, std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < m; i++) {
        uint64_t o = target_vectors[i];
        if (o == 0)
            continue;
        for (size_t j = 0; j < full_basis.size(); j++)
            if (full_basis[j] == o) {
                result.method.outputs[i] = j;
                break;
            }

        assert(result.method.outputs[i] != std::numeric_limits<size_t>::max());
    }

    assert(B.empty());
    result.additions_before = gf2::gp::naive_additions(G.matrix, m);
    result.additions_after = result.method.additions.size();

    return result;
}

} // namespace slp::mip
#else

#include <stdexcept>

namespace slp::mip {

Result run_MIP(const Z2Matrix &, const size_t, const size_t) {
    throw std::runtime_error("MIP backend unavailable: rebuild with "
                             "USE_MATHOPT=1 and OR-Tools installed");
}

} // namespace slp::mip

#endif
