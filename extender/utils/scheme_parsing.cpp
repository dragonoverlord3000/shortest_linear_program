#include <bits/stdc++.h>
using namespace std;
namespace fs = std::filesystem;

// collects all the scheme-files
vector<string> collect_recursive(const fs::path &root_dir) {
    vector<string> out;
    for (const auto &e : fs::recursive_directory_iterator(root_dir)) {
        if (e.is_regular_file())
            out.push_back(e.path().string());
    }
    return out;
}

// get all lines from a given file
static inline vector<string> read_all_lines(const fs::path &file) {
    ifstream in(file);
    if (!in)
        throw runtime_error("Failed to open: " + file.string());

    vector<string> lines;
    string line;
    while (getline(in, line))
        lines.push_back(line);
    return lines;
}

// small helper to split line
static inline vector<string> split_char(const string &s, char delim) {
    vector<string> parts;
    string cur;
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
static inline string trim(const string &s) {
    size_t a = 0, b = s.size();
    while (a < b && isspace(static_cast<unsigned char>(s[a])))
        a++;
    while (b > a && isspace(static_cast<unsigned char>(s[b - 1])))
        b--;
    return s.substr(a, b - a);
}

// parse a row into a vector<int>
static inline vector<int> parse_int_row_space_separated(const std::string &s) {
    vector<int> row;
    std::istringstream iss(s);
    int x;
    while (iss >> x)
        row.push_back(x); // skips multiple spaces automatically
    return row;
}
