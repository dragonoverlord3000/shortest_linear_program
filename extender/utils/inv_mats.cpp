#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

// finds the GF2 inverse of the matrix B
vector<vector<int>> inv_mat_mod2(vector<vector<int>> B) {
    int n = B.size();
    vector<vector<int>> B_inv(n, vector<int>(n));
    for (int i = 0; i < n; i++)
        B_inv[i][i] = 1;

    for (int row = 0; row < n; row++) {
        int pivot_row = -1;
        for (int possible_pivot = row; possible_pivot < n; possible_pivot++) {
            if (B[possible_pivot][row]) {
                pivot_row = possible_pivot;
                break;
            }
        }

        // assumes B invertible, i.e. pivot_row = -1 should never be possible
        // here
        assert(pivot_row != -1);
        swap(B[pivot_row], B[row]);
        swap(B_inv[pivot_row], B_inv[row]);

        for (int r = 0; r < n; r++) {
            if (r != row and B[r][row]) {
                for (int k = 0; k < n; k++) {
                    B[r][k] ^= B[row][k];
                    B_inv[r][k] ^= B_inv[row][k];
                }
            }
        }
    }

    return B_inv;
}

// a matrix that is invertible over the integers is an unimodular matrix and is
// invertible over Z2, if we assume that our scheme comes from such a matrix,
// then we can try all {-1,0,1} unimodular extensions of our GF2 matrix as a
// heuristic. https://en.wikipedia.org/wiki/Unimodular_matrix will return a
// fraction (numerator / denominator)
class Fraction {
  public:
    long long num;
    long long den;

    Fraction() : num(0), den(1){};

    Fraction(long long num, long long den) : num(num), den(den) { reduce(); }

    void reduce() {
        if (num) {
            long long g = gcd(llabs(num), llabs(den));
            num /= g;
            den /= g;
        } else
            den = 1;
        if (den < 0) {
            den *= -1;
            num *= -1;
        }
    }

    void operator*=(long long f) {
        num *= f;
        reduce();
    }

    void operator/=(Fraction other) {
        num *= other.den;
        den *= other.num;
        reduce();
    }

    friend Fraction operator+(Fraction a, Fraction b) {
        long long numerator = a.num * b.den + b.num * a.den;
        long long denominator = a.den * b.den;
        return Fraction(numerator, denominator);
    }

    friend Fraction operator-(Fraction a, Fraction b) {
        long long numerator = a.num * b.den - b.num * a.den;
        long long denominator = a.den * b.den;
        return Fraction(numerator, denominator);
    }

    friend Fraction operator*(Fraction a, Fraction b) {
        long long numerator = a.num * b.num;
        long long denominator = a.den * b.den;
        return Fraction(numerator, denominator);
    }

    friend Fraction operator/(Fraction a, Fraction b) {
        long long numerator = a.num * b.den;
        long long denominator = a.den * b.num;
        return Fraction(numerator, denominator);
    }
};

Fraction compute_det(const vector<vector<int>> &_mat) {
    int n = _mat.size();
    Fraction det(1, 1);
    vector<vector<Fraction>> mat(n, vector<Fraction>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mat[i][j] = {_mat[i][j], 1};

    // gaussian elimination
    for (int row = 0; row < n; row++) {
        int pivot_row = -1;
        for (int r = row; r < n; r++) {
            if (mat[r][row].num) {
                pivot_row = r;
                break;
            }
        }
        if (pivot_row == -1)
            return Fraction(0, 1);
        swap(mat[pivot_row], mat[row]);
        if (pivot_row != row)
            det *= -1;

        Fraction pivot = mat[row][row];

        for (int r = row + 1; r < n; r++) {
            Fraction mult_factor = mat[r][row] / pivot;
            for (int k = row; k < n; k++)
                mat[r][k] = mat[r][k] - mult_factor * mat[row][k];
        }
    }

    // computing the determinant
    for (int i = 0; i < n; i++)
        det = det * mat[i][i];

    return det;
}

// compute the integer inverse of a proven unimodular matrix
vector<vector<int>> int_inv(vector<vector<int>> &mat) {
    int n = mat.size();
    vector<vector<Fraction>> mat_frac(n, vector<Fraction>(n));
    vector<vector<Fraction>> mat_inv_frac(n, vector<Fraction>(n));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mat_frac[i][j] = {mat[i][j], 1};

    for (int i = 0; i < n; i++)
        mat_inv_frac[i][i] = Fraction(1, 1);

    // compute inverse
    for (int row = 0; row < n; row++) {
        int pivot_row = -1;
        for (int r = row; r < n; r++) {
            if (mat_frac[r][row].num) {
                pivot_row = r;
                break;
            }
        }
        // should never happen, just a sanity check
        assert(pivot_row != -1);
        swap(mat_frac[row], mat_frac[pivot_row]);
        swap(mat_inv_frac[row], mat_inv_frac[pivot_row]);

        Fraction d = mat_frac[row][row];
        for (int k = 0; k < n; k++) {
            mat_frac[row][k] = mat_frac[row][k] / d;
            mat_inv_frac[row][k] = mat_inv_frac[row][k] / d;
        }

        for (int r = 0; r < n; r++) {
            if (r == row)
                continue;
            Fraction mult_factor = mat_frac[r][row];
            for (int k = 0; k < n; k++) {
                mat_frac[r][k] =
                    mat_frac[r][k] - mult_factor * mat_frac[row][k];
                mat_inv_frac[r][k] =
                    mat_inv_frac[r][k] - mult_factor * mat_inv_frac[row][k];
            }
        }
    }

    // port to int (unimodular implies this can be done)
    vector<vector<int>> mat_inv(n, vector<int>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            // sanity check
            assert((llabs(mat_inv_frac[i][j].num) %
                    llabs(mat_inv_frac[i][j].den)) == 0);
            mat_inv[i][j] = mat_inv_frac[i][j].num / mat_inv_frac[i][j].den;
        }

    return mat_inv;
}

// will find all unimodular {-1,0,1} matrices that reduce to the GF2 mat found
// note that this method relies on the fact that bit_p in our heuristic is low
// we return these matrices
vector<vector<vector<int>>> unimods(vector<vector<int>> &mat,
                                    bool verbose = false) {
    int n = mat.size();
    vector<pair<int, int>> set_cells;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (mat[i][j])
                set_cells.push_back({i, j});

    if (verbose)
        cout << "number of set cells in basis change matrix: "
             << set_cells.size() << endl;

    // note that it will fail if set_cells is greater than 63, but this should
    // never happen as per the bit_p being reasonably low
    uint64_t one = 1;
    uint64_t num_attempts = one << set_cells.size();

    // list of unimodular matrices that reduce to our GF2 matrix
    vector<vector<vector<int>>> ans;

    for (uint64_t bitmask = 0; bitmask < num_attempts; bitmask++) {
        vector<vector<int>> B(n, vector<int>(n));
        for (size_t shift = 0; shift < set_cells.size(); shift++) {
            int i = set_cells[shift].first;
            int j = set_cells[shift].second;
            B[i][j] = bitmask & (one << shift) ? 1 : -1;
        }

        Fraction det = compute_det(B);
        det.reduce();
        if (llabs(det.num) != 1 || llabs(det.den) != 1)
            continue;
        ans.push_back(B);
    }

    return ans;
}

void print_mat(vector<vector<int>> &mat) {
    int n = mat.size();
    cout << "[" << endl;
    for (int i = 0; i < n; i++) {
        cout << "[";
        for (int j = 0; j < n; j++)
            cout << mat[i][j] << ", "[j == n - 1];
        cout << "]" << endl;
    }
    cout << "]" << endl;
}

/*
// finds the general real number inverse of B
vector<vector<double>> inv_mat_general(vector<vector<int>> B) {
    // TODO
    return {};
}

 * This is just for some simple testing
vector<vector<int>> matmul(vector<vector<int>>& A, vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n));
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            for(int k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

int main() {
    vector<vector<int>> test_mat = {
        {1,0,1,0},
        {1,1,0,0},
        {0,0,1,0},
        {0,0,1,1}
    };
    cout << "input matrix Z2: " << endl;
    print_mat(test_mat);

    vector<vector<vector<int>>> unimod_extension_mats = unimod_invs(test_mat);
    cout << "found " << unimod_extension_mats.size() << " matrices" << endl;
    if (!unimod_extension_mats.empty()) {
        cout << "example inv: " << endl;
        print_mat(unimod_extension_mats[0]);
        vector<vector<int>> inv0 = int_inv(unimod_extension_mats[0]);
        cout << "example og: " << endl;
        print_mat(inv0);
        vector<vector<int>> product = matmul(unimod_extension_mats[0], inv0);
        cout << "example general real number product: " << endl;
        print_mat(product);
    }

    return 0;
}
*/
