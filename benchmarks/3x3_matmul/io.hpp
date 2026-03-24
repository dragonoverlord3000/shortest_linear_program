#pragma once

#include <filesystem>
#include <vector>

namespace io {
namespace fs = std::filesystem;

enum class SchemeType { W, U, V };

std::vector<std::string> collect_recursive(const fs::path &root_dir);
std::vector<std::string> read_all_lines(const fs::path &file);
std::vector<std::string> split_char(const std::string &s, char delim);
std::string trim(const std::string &s);
std::vector<uint64_t> parse_int_row_space_separated(const std::string &s);
void update_G_W(std::vector<uint64_t> &G, const std::string &s_row, int i, int col_idx);
void update_G_UV(std::vector<uint64_t> &G, const std::string &s_row, int i, int row_idx);
std::vector<uint64_t> parse_one_file(const fs::path &file, const std::string &type);

} // namespace io
