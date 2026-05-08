#include "io.hpp"

#include <cassert>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace io {
// READ ALL SCHEMES START
// collects all the scheme-files
std::vector<std::string> collect_recursive(const fs::path &root_dir) {
    std::vector<std::string> out;
    for (const auto &e : fs::recursive_directory_iterator(root_dir)) {
        if (e.is_regular_file())
            out.push_back(e.path().string());
    }
    return out;
}

// get all lines from a given file
std::vector<std::string> read_all_lines(const fs::path &file) {
    std::ifstream in(file);
    if (!in)
        throw std::runtime_error("Failed to open: " + file.string());

    std::vector<std::string> lines;
    std::string line;
    while (getline(in, line))
        lines.push_back(line);
    return lines;
}

// small helper to split line
std::vector<std::string> split_char(const std::string &s, char delim) {
    std::vector<std::string> parts;
    std::string cur;
    for (char c : s) {
        if (c == delim) {
            parts.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    parts.push_back(cur);
    return parts;
}

// remove prefix and suffix whitespace
std::string trim(const std::string &s) {
    size_t a = 0, b = s.size();
    while (a < b && isspace(static_cast<unsigned char>(s[a])))
        a++;
    while (b > a && isspace(static_cast<unsigned char>(s[b - 1])))
        b--;
    return s.substr(a, b - a);
}

// parse a row into a vector<int>
std::vector<uint64_t> parse_int_row_space_separated(const std::string &s) {
    std::vector<uint64_t> row;
    std::istringstream iss(s);
    int x;
    while (iss >> x)
        row.push_back(x); // skips multiple spaces automatically
    return row;
}
// READ ALL SCHEMES END

// PARSE ONE FILE START
void update_G_W(std::vector<uint64_t> &G, const std::string &s_row, int i,
                int col_idx) {
    std::vector<uint64_t> row = parse_int_row_space_separated(s_row);
    for (int j = 0; j < 3; j++) {
        int row_idx = 3 * i + j;
        if (row[j])
            G[col_idx] |= 1u << row_idx;
    }
}

void update_G_UV(std::vector<uint64_t> &G, const std::string &s_row, int i,
                 int row_idx) {
    std::vector<uint64_t> row = parse_int_row_space_separated(s_row);
    for (int j = 0; j < 3; j++) {
        int col_idx = 3 * i + j;
        if (row[j])
            G[col_idx] |= 1u << row_idx;
    }
}

std::vector<uint64_t> parse_one_file_3x3_matmul(const fs::path &file,
                                                const std::string &type) {
    auto raw = read_all_lines(file);

    std::vector<uint64_t> G(type == "W" ? 23 : 9, 0);
    size_t main_idx = 0;

    for (int col_idx = 0; col_idx < 23; ++col_idx) {
        for (int i = 0; i < 3; i++) {
            const std::string &row = raw[main_idx];
            auto split_row = split_char(row, '|');
            if (type == "W") {
                update_G_W(G, trim(split_row[2]), i, col_idx);
            } else {
                update_G_UV(G, trim(split_row[type == "V"]), i, col_idx);
            }
            main_idx++;
        }
        // skip separator
        main_idx++;
    }
    return G;
}

std::tuple<size_t, size_t, std::vector<uint64_t>>
parse_one_file_crypt(const fs::path &file) {
    auto raw = read_all_lines(file);

    std::string row1 = trim(raw[1]);
    std::vector<uint64_t> m_n = parse_int_row_space_separated(row1);
    size_t m = m_n[0], n = m_n[1];

    std::vector<uint64_t> G(n, 0);
    for (size_t i = 0; i < m; i++) {
        std::string row = trim(raw[i + 2]);
        std::vector<uint64_t> r = parse_int_row_space_separated(row);
        assert(r.size() == n);
        for (size_t j = 0; j < n; j++)
            G[j] |= r[j] << i;
    }

    return {m, n, G};
}

std::tuple<size_t, size_t, std::vector<uint64_t>>
parse_one_file_struct(const fs::path &file) {
    auto raw = read_all_lines(file);

    std::string row0 = trim(raw[0]);
    std::vector<uint64_t> m_n = parse_int_row_space_separated(row0);
    size_t m = m_n[0], n = m_n[1];

    std::vector<uint64_t> G(n, 0);
    for (size_t i = 0; i < m; i++) {
        std::string row = trim(raw[i + 1]);
        std::vector<uint64_t> r = parse_int_row_space_separated(row);
        assert(r.size() == n);
        for (size_t j = 0; j < n; j++)
            G[j] |= r[j] << i;
    }

    return {m, n, G};
}
// PARSE ONE FILE END
} // namespace io
