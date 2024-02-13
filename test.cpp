/*
    OUTPUT:
        C matrix solution from function naive_row_matmul:
        Time taken by function: 441200 microseconds
        C solution from function avx256_row_matmul:
        Time taken by function: 225400 microseconds

    Where,

    constexpr int m = 70;
    constexpr int n = 90;
    constexpr int k = 40;

    Solving matrix-matrix multiple: 
    C[m][k] = C[m][k] + A[m][n] * B[n][k]
*/

// #include <lapack.h>
// #include "lapack.h"
#include "matmul.hpp"
#include <chrono>

template <typename T, std::size_t Rows, std::size_t Cols>
using Matrix = matrix::Matrix<T,Rows,Cols>;

constexpr int m = 70;
constexpr int n = 90;
constexpr int k = 40;

// Needed to find leading dimensional value at compile time
// Only useful for padding the arrays
// Would lead to some minor optimizations over handling fringe cases
#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

constexpr int mbuffer = (m % 4); 
constexpr int nbuffer = (n % 4); 
constexpr int kbuffer = (k % 4); 
constexpr int maxbuffer = max(mbuffer, max(nbuffer, kbuffer)); 
constexpr int ldmax = max(m + mbuffer, max(n + nbuffer, k + kbuffer));

int main() {
    //std::uniform_real_distribution<> dis = matrix::init_rand_gen(gen,0.0, 10.0);

    Matrix<double,m,n> A; 
    matrix::fill_matrix<double,m,n>(1.0,A);

    Matrix<double,n,k> B;
    matrix::fill_matrix<double,n,k>(1.0,B);

    Matrix<double,m,k> C;
    matrix::fill_matrix<double,m,k>(0.0,C);

    // matrix::print_matrix<double,m,n>(A);
    // std::cout << std::endl;

    // matrix::print_matrix<double,n,k>(B);
    // std::cout << std::endl;

    // matrix::print_matrix<double,m,k>(C);
    // std::cout << std::endl;

    /* IGNORE BUFFERING MATRICES FOR NOW */
    // Matrix<double,ldmax,ldmax> C_buf = matrix::create_pad_matrix<double,m,k,ldmax>(C);
    // matrix::print_matrix<double,ldmax,ldmax>(C_buf);

    // std::cout << "ldmax = " << ldmax << std::endl;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    naive_row_matmul(m,n,k,A,B,C);
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = endTime - startTime;

    std::cout << "C matrix solution from function naive_row_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
    // matrix::print_matrix<double,m,k>(C);

    matrix::fill_matrix<double,m,k>(0.0,C);

    startTime = std::chrono::high_resolution_clock::now();
    avx256_row_matmul(m,n,k,A,B,C);
    endTime = std::chrono::high_resolution_clock::now();
    duration = endTime - startTime;

    std::cout << "C solution from function avx256_row_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
    // matrix::print_matrix<double,m,k>(C);
}
