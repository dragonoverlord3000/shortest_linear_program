#include "slp/types.hpp"

#include <chrono>
#include <iostream>
#include <vector>

namespace slp::gf2 {
void validate_method(const AdditionMethod &method, size_t n, const char *where);

void toposorter(AdditionMethod &addition_method, const size_t n);

// small duration utility
template <class Duration>
double
remaining_seconds(const std::chrono::time_point<std::chrono::steady_clock,
                                                Duration> &deadline) {
    const auto now = std::chrono::steady_clock::now();
    return std::max(0.0, std::chrono::duration<double>(deadline - now).count());
}

// small print helper
template <typename T> inline void print_vec(std::vector<T> &arr) {
    std::cout << "{";
    for (size_t i = 1; i < arr.size(); i++)
        std::cout << arr[i - 1] << ", ";

    if (!arr.empty())
        std::cout << arr.back();
    std::cout << "}" << std::endl;
}
} // namespace slp::gf2
