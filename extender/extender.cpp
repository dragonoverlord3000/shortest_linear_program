#include <algorithm>
#include <cassert>
#include <cctype>
#include <ctime> // clock, CLOCKS_PER_SEC
#include <filesystem>
#include <fstream>
#include <iostream> // cout, cin, ios_base
#include <map>
#include <set>
#include <stdexcept> // runtime_error
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// External
#include "nlohmann/json.hpp"
using json = nlohmann::json;

// Internal imports
#include "utils/inv_mats.cpp"
#include "utils/random.cpp"
#include "utils/scheme_parsing.cpp"

using namespace std;

namespace fs = std::filesystem;

template <typename T> void print_v(T iterable) {
    cout << "{";
    for (auto v : iterable)
        cout << v << ",";
    cout << "\b}";
}

// START OF READING GIVEN SCHEME
void update_G_W(vector<vector<int>> &G, const string &s_row, int i,
                int col_idx) {
    vector<int> row = parse_int_row_space_separated(s_row);
    if (row.size() != 3)
        throw runtime_error("Expected 3 ints in row: '" + s_row + "'");

    for (int j = 0; j < 3; j++) {
        int row_idx = 3 * i + j;
        G[row_idx][col_idx] = row[j];
    }
}

void update_G_UV(vector<vector<int>> &G, const string &s_row, int i,
                 int row_idx) {
    vector<int> row = parse_int_row_space_separated(s_row);
    if (row.size() != 3)
        throw runtime_error("Expected 3 ints in row: '" + s_row + "'");

    for (int j = 0; j < 3; j++) {
        int col_idx = 3 * i + j;
        G[row_idx][col_idx] = row[j];
    }
}

// Parse one file and get out U,V or W
vector<vector<int>> parse_one_file(const fs::path &file, const string &type) {
    auto raw = read_all_lines(file);

    int m = (type == "W") ? 9 : 23;
    int n = (type == "W") ? 23 : 9;
    vector<vector<int>> G(m, vector<int>(n));
    size_t main_idx = 0;

    for (int col_idx = 0; col_idx < 23; ++col_idx) {
        for (int i = 0; i < 3; i++) {
            if (main_idx >= raw.size())
                throw runtime_error("Parse OOB in " + file.string() +
                                    " at line " + to_string(main_idx));

            const string &row = raw[main_idx];
            auto split_row = split_char(row, '|');
            if (split_row.size() < 3)
                throw runtime_error("Bad row (missing |) in " + file.string() +
                                    " at line " + to_string(main_idx) + ": " +
                                    row);

            int type_idx = 0;
            if (type == "V")
                type_idx = 1;
            if (type == "W")
                type_idx = 2;

            if (type == "W")
                update_G_W(G, trim(split_row[type_idx]), i, col_idx);
            else
                update_G_UV(G, trim(split_row[type_idx]), i, col_idx);
            ++main_idx;
        }
        // skip separator if present
        if (main_idx < raw.size())
            ++main_idx;
    }
    return G;
}

// G <- BG
inline void _change_basis_W(const vector<vector<int>> &G,
                            const vector<vector<int>> &B,
                            vector<vector<int>> &out) {
    size_t m = G.size();
    size_t n = G[0].size();
    out.assign(m, vector<int>(n, 0));

    assert(m == B[0].size());

    for (size_t i = 0; i < B.size(); i++)
        for (size_t j = 0; j < n; j++)
            for (size_t k = 0; k < m; k++)
                out[i][j] += B[i][k] * G[k][j];
}

// G <- GB
inline void _change_basis_UV(const vector<vector<int>> &G,
                             const vector<vector<int>> &B,
                             vector<vector<int>> &out) {
    size_t m = G.size();
    size_t n = G[0].size();
    out.assign(m, vector<int>(n, 0));

    assert(n == B.size());

    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < B[0].size(); j++)
            for (size_t k = 0; k < n; k++)
                out[i][j] += G[i][k] * B[k][j];
}
// END OF READING GIVEN SCHEME

// START OF READING FOUND ADDITION SCHEME
class Node {
  public:
    string varname;
    vector<int> vec;
    Node(string &varname, int n) : varname(varname) {
        vec.assign(n, 0);

        // the building blocks are the singleton elements, all other elements
        // should be built from these
        if (varname[0] == 's') {
            // note that all vars are of the form x_i for x in s,u,v,w and i an
            // integer the first n elements are singletons
            int idx = stoi(varname.substr(2)) - 1;
            if (idx < n)
                vec[idx] = 1;
        }
    }
};

// small helper
vector<string> splitter(string &s) {
    vector<string> out;
    string cur = "";
    for (char c : s) {
        if (c == ' ') {
            if (!cur.empty())
                out.push_back(cur);
            cur.clear();
        } else
            cur.push_back(c);
    }
    if (!cur.empty())
        out.push_back(cur);
    return out;
}
static string normalized_path_string(string s) {
    replace(s.begin(), s.end(), '\\', '/');
    return s;
}

static bool same_scheme_path(const string &a_raw, const string &b_raw) {
    string a = normalized_path_string(a_raw);
    string b = normalized_path_string(b_raw);

    return a == b || a.ends_with(b) || b.ends_with(a);
}

static int mod2(int x) {
    x %= 2;
    if (x < 0)
        x += 2;
    return x;
}

static string vec_key(const vector<int> &v) {
    string s;
    s.reserve(v.size());
    for (int x : v)
        s.push_back(mod2(x) ? '1' : '0');
    return s;
}

static vector<vector<int>> bitmasks_to_matrix(const json &masks) {
    vector<vector<int>> B;
    int n = static_cast<int>(masks.size());
    B.assign(n, vector<int>(n, 0));

    for (int i = 0; i < n; ++i) {
        uint64_t mask = masks[i].get<uint64_t>();
        for (int j = 0; j < n; ++j)
            B[i][j] = static_cast<int>((mask >> j) & 1ULL);
    }

    return B;
}

static vector<string> make_baselines(int n) {
    vector<string> out;
    for (int i = 1; i <= n; ++i)
        out.push_back("s_" + to_string(i));
    return out;
}

static vector<vector<int>> simulate_best_method_vectors(const json &method,
                                                        int n) {
    vector<vector<int>> vecs;

    for (int i = 0; i < n; ++i) {
        vector<int> e(n, 0);
        e[i] = 1;
        vecs.push_back(e);
    }

    for (const auto &op : method) {
        int a = op[0].get<int>();
        int b = op[1].get<int>();

        if (a < 0 || b < 0 || a >= static_cast<int>(vecs.size()) ||
            b >= static_cast<int>(vecs.size())) {
            throw runtime_error("Bad best_method index");
        }

        vector<int> v(n, 0);
        for (int i = 0; i < n; ++i)
            v[i] = vecs[a][i] ^ vecs[b][i];

        vecs.push_back(v);
    }

    return vecs;
}

static vector<string> best_method_to_intermediates(const json &method, int n) {
    vector<string> out;

    for (int i = 0; i < static_cast<int>(method.size()); ++i) {
        int dst = n + i + 1;                 // old node names are 1-based
        int a = method[i][0].get<int>() + 1; // new indices are 0-based
        int b = method[i][1].get<int>() + 1;

        out.push_back("s_" + to_string(dst) + " s_" + to_string(a) + " s_" +
                      to_string(b));
    }

    return out;
}

static vector<vector<int>>
apply_basis_change_to_G(const vector<vector<int>> &G0,
                        const vector<vector<int>> &basis_change,
                        const string &type) {
    vector<vector<int>> G;

    if (type == "W")
        _change_basis_W(G0, basis_change, G);
    else
        _change_basis_UV(G0, basis_change, G);

    for (auto &row : G)
        for (int &x : row)
            x = mod2(x);

    return G;
}

static vector<string> infer_outputs_from_best_method(
    const json &method, const vector<vector<int>> &G0,
    const vector<vector<int>> &basis_change, const string &type, int n) {
    vector<vector<int>> vecs = simulate_best_method_vectors(method, n);

    unordered_map<string, int> key2idx;
    for (int i = 0; i < static_cast<int>(vecs.size()); ++i)
        key2idx.emplace(vec_key(vecs[i]), i);

    vector<vector<int>> G = apply_basis_change_to_G(G0, basis_change, type);

    string prefix = type;
    transform(prefix.begin(), prefix.end(), prefix.begin(),
              [](unsigned char c) { return tolower(c); });

    vector<string> outputs;

    for (int i = 0; i < static_cast<int>(G.size()); ++i) {
        string key = vec_key(G[i]);

        if (!key2idx.count(key)) {
            throw runtime_error("Could not infer output node for " + type +
                                "_" + to_string(i + 1));
        }

        int node_idx = key2idx[key] + 1; // old names are 1-based
        outputs.push_back(prefix + "_" + to_string(i + 1) + " s_" +
                          to_string(node_idx));
    }

    return outputs;
}

// outputs basis change matrix (Z2) and edge list
pair<vector<vector<int>>, vector<pair<Node *, Node *>>>
read_add_scheme(string &addition_scheme_file_path, string &scheme, string &type,
                const vector<vector<int>> &G0) {
    ifstream f(addition_scheme_file_path);
    json addition_schemes = json::parse(f);

    bool found = false;
    vector<vector<int>> basis_change;
    vector<string> baselines, intermediates, outputs;

    // benchmark format - designed specifically for rank-23 3x3 matmul
    for (const auto &benchmark_result : addition_schemes["results"]) {
        const auto &details = benchmark_result["details"];

        if (!details.contains("schemes"))
            continue;

        for (const auto &scheme_entry : details["schemes"]) {
            string json_scheme = scheme_entry["scheme_name"].get<string>();

            if (!same_scheme_path(json_scheme, scheme))
                continue;

            const auto &results = scheme_entry["results"];

            if (!results.contains(type))
                continue;

            const auto &r = results[type];

            found = true;

            basis_change = bitmasks_to_matrix(r["best_basis_change"]);

            int n = (type == "W") ? 23 : 9;

            baselines = make_baselines(n);
            intermediates = best_method_to_intermediates(r["best_method"], n);
            outputs = infer_outputs_from_best_method(r["best_method"], G0,
                                                     basis_change, type, n);

            break;
        }

        if (found)
            break;
    }

    if (!found) {
        throw runtime_error("Could not find scheme/type in addition JSON: " +
                            scheme + " / " + type);
    }

    // the length of the row vectors
    int n = (type == "W") ? 23 : 9;

    unordered_map<string, Node *> nodename2node;
    vector<pair<Node *, Node *>> edges;
    vector<vector<string>> t = {baselines, intermediates, outputs};

    for (vector<string> &ss : t) {
        for (string &s : ss) {
            vector<string> elements = splitter(s);
            string to_build = elements[0];

            if (!nodename2node.count(to_build))
                nodename2node[to_build] = new Node(to_build, n);

            for (size_t i = 1; i < elements.size(); i++) {
                string atom = elements[i];

                if (!nodename2node.count(atom))
                    nodename2node[atom] = new Node(atom, n);

                edges.push_back({nodename2node[atom], nodename2node[to_build]});
            }
        }
    }

    return {basis_change, edges};
}
// END OF READING FOUND ADDITION SCHEME

// START OF BACKTRACKING SOLVER
pair<bool, vector<int>>
dfs_get_row(Node *node, unordered_map<Node *, vector<pair<Node *, int>>> &adj) {
    if (adj[node].empty()) {
        // a small sanity check
        size_t idx = stoi(node->varname.substr(2)) - 1;
        assert(idx < node->vec.size());
        return {false, node->vec};
    }

    vector<int> out(node->vec.size(), 0);
    bool undecided = false;

    for (auto &[neigh, f] : adj[node]) {
        undecided |= f == 0;
        if (undecided)
            break;
        auto [sub_undecided, sub_vec] = dfs_get_row(neigh, adj);
        undecided |= sub_undecided;
        for (size_t i = 0; i < out.size(); i++)
            out[i] += f * sub_vec[i];
    }

    return {undecided, out};
}
// returns true if any of the edges related to row_node is undetermined (0
// weight) or it conicides with expectation as per G, otherwise false
bool row_node_checks_out(Node *row_node,
                         unordered_map<Node *, vector<pair<Node *, int>>> &adj,
                         vector<vector<int>> &G) {
    auto [undetermined, vec_ans] = dfs_get_row(row_node, adj);
    if (undetermined)
        return true;
    int idx = stoi(row_node->varname.substr(2)) - 1;
    for (size_t i = 0; i < vec_ans.size(); i++)
        if (G[idx][i] != vec_ans[i])
            return false;
    return true;
}

// returns false if any u,v,w all of whos edges have defined factors does not
// give result matching G
bool is_valid_so_far(vector<pair<Node *, Node *>> &edges, vector<int> &factors,
                     vector<vector<int>> &G, Node *row_node) {
    // first element is neighborg node, second is factor, or 0 if no factor for
    // this edge yet note that we swap the order here so that it goes from final
    // (u,v,w) to lowest building blovk (s_1, s_2, s_3, s_4, ...)
    unordered_map<Node *, vector<pair<Node *, int>>> adj;
    for (size_t i = 0; i < factors.size(); i++) {
        auto [node_from, node_to] = edges[i];
        int factor = factors[i];
        adj[node_to].push_back({node_from, factor});
    }
    for (size_t i = factors.size(); i < edges.size(); i++) {
        auto [node_from, node_to] = edges[i];
        adj[node_to].push_back({node_from, 0});
    }

    if (!row_node_checks_out(row_node, adj, G))
        return false;

    return true;
}

bool bt(size_t at, vector<pair<Node *, Node *>> &edges, vector<int> &factors,
        vector<vector<int>> &G, unordered_map<int, Node *> &idx2row_node) {
    if (at == edges.size())
        return true;
    vector<int> possible = {-1, 1};
    for (int v : possible) {
        factors.push_back(v);
        if (idx2row_node.count(at) &&
            !is_valid_so_far(edges, factors, G, idx2row_node[at])) {
            factors.pop_back();
            continue;
        }
        if (bt(at + 1, edges, factors, G, idx2row_node))
            return true;
        factors.pop_back();
    }

    return false;
}

//
int dfs_get_max_idx(Node *node, unordered_map<Node *, vector<Node *>> &adj,
                    map<pair<Node *, Node *>, int> &edge2idx) {
    int max_idx = -1;
    for (Node *neigh : adj[node]) {
        pair<Node *, Node *> edge = {neigh, node};
        max_idx = max(max_idx, edge2idx[edge]);
        max_idx = max(max_idx, dfs_get_max_idx(neigh, adj, edge2idx));
    }

    return max_idx;
}

void dfs_fill_sorted_edges(Node *node,
                           unordered_map<Node *, vector<Node *>> &adj,
                           vector<pair<Node *, Node *>> &sorted_edges,
                           set<pair<Node *, Node *>> &seen_edges) {
    for (Node *neigh : adj[node]) {
        pair<Node *, Node *> edge = {neigh, node};
        if (!seen_edges.count(edge)) {
            sorted_edges.push_back(edge);
            seen_edges.insert(edge);
        }
        dfs_fill_sorted_edges(neigh, adj, sorted_edges, seen_edges);
    }
}

vector<pair<Node *, Node *>> edge_sorter(vector<pair<Node *, Node *>> &edges) {
    // first find row_nodes and initial adjacency list (not it goes from top
    // (u,v,w) down)
    unordered_set<Node *> row_nodes;
    unordered_map<Node *, vector<Node *>> _adj;
    for (auto [node_from, node_to] : edges) {
        _adj[node_to].push_back(node_from);
        if (node_from->varname[0] != 's')
            row_nodes.insert(node_from);
        if (node_to->varname[0] != 's')
            row_nodes.insert(node_to);
    }

    // sort the edges to waste as little work as possible
    vector<pair<Node *, Node *>> sorted_edges;
    set<pair<Node *, Node *>> seen_edges;
    for (Node *node : row_nodes)
        dfs_fill_sorted_edges(node, _adj, sorted_edges, seen_edges);

    // sanity check
    assert(sorted_edges.size() == edges.size());

    return sorted_edges;
}

unordered_map<int, Node *>
get_idx2row_node(const vector<pair<Node *, Node *>> &edges) {
    unordered_set<Node *> row_nodes;
    map<pair<Node *, Node *>, int> edge2idx;
    // note that we build adj to go the opposite way of edges i.e. from top
    // (u,v,w) to bottom (s_1,s_2,s_3,...)
    unordered_map<Node *, vector<Node *>> adj;
    for (size_t i = 0; i < edges.size(); i++) {
        auto [node_from, node_to] = edges[i];
        adj[node_to].push_back(node_from);
        edge2idx[{node_from, node_to}] = i;
        if (node_from->varname[0] != 's')
            row_nodes.insert(node_from);
        if (node_to->varname[0] != 's')
            row_nodes.insert(node_to);
    }

    unordered_map<int, Node *> idx2row_node;
    for (Node *row_node : row_nodes) {
        int idx = dfs_get_max_idx(row_node, adj, edge2idx);
        assert(idx >= 0);
        idx2row_node[idx] = row_node;
    }

    return idx2row_node;
}

pair<bool, vector<tuple<string, string, int>>>
backtracker(vector<vector<int>> &G, vector<pair<Node *, Node *>> &edges,
            unordered_map<int, Node *> &idx2row_node) {
    vector<int> factors;
    bool works = bt(0, edges, factors, G, idx2row_node);

    vector<tuple<string, string, int>> solution;
    for (size_t i = 0; i < factors.size(); i++)
        solution.push_back(
            {edges[i].first->varname, edges[i].second->varname, factors[i]});

    return {works, solution};
}
// END OF BACKTRACKING SOLVER

void solve(string &scheme, string &type, string &addition_scheme_file_path) {
    string file = scheme;
    vector<vector<int>> G0 = parse_one_file(file, type);
    auto [basis_change, edges] =
        read_add_scheme(addition_scheme_file_path, scheme, type, G0);
    edges = edge_sorter(edges);
    unordered_map<int, Node *> idx2row_node = get_idx2row_node(edges);

    // find all unimodular {-1, 0, 1} extensions
    vector<vector<vector<int>>> all_possible_basis_change =
        unimods(basis_change);

    cout << "num unimodular extensions found: "
         << all_possible_basis_change.size() << endl;

    int cnt = 0;
    for (vector<vector<int>> &B : all_possible_basis_change) {
        cnt++;
        cout << "cnt: " << cnt << endl;
        vector<vector<int>> G;
        if (type == "W")
            _change_basis_W(G0, B, G);
        else
            _change_basis_UV(G0, B, G);

        // call backtrack solver and log result
        auto [works, solution] = backtracker(G, edges, idx2row_node);
        if (!works) {
            cout << "no solution found {-1, 0, 1}" << endl;
            continue;
        }
        cout << "solution found" << endl;
        cout << "basis change matrix: " << endl;
        print_mat(B);
        cout << "factors: " << endl;
        for (auto [node_from_name, node_from_to, factor] : solution)
            cout << node_from_name << " " << node_from_to << " " << factor
                 << endl;
        return;
    }
}

int main(int argc, char **argv) {
    string type;
    string scheme;
    string addition_scheme_file_path;

    if (argc == 4) {
        type = argv[1];
        scheme = argv[2];
        addition_scheme_file_path = argv[3];
    } else {
        cin >> type;
        cin >> scheme;
        cin >> addition_scheme_file_path;
    }

    solve(scheme, type, addition_scheme_file_path);

    return 0;
}
