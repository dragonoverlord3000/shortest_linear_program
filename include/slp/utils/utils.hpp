#include <iostream>
#include <vector>

template <typename T> inline void print_vec(std::vector<T> &arr) {
    std::cout << "{";
    for (size_t i = 1; i < arr.size(); i++)
        std::cout << arr[i - 1] << ", ";

    if (!arr.empty())
        std::cout << arr.back();
    std::cout << "}" << std::endl;
}
