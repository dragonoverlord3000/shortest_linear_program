#include <iostream>
#include <vector>

void print_mat(std::vector<std::vector<int>> &mat) {
    int n = mat.size();
    std::cout << "[" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "[";
        for (int j = 0; j < n; j++)
            std::cout << mat[i][j] << ", "[j == n - 1];
        std::cout << "]" << std::endl;
    }
    std::cout << "]" << std::endl;
}
