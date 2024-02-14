/*
    std::array OUTPUT:

        C matrix solution from function naive_row_matmul:
        Time taken by function: 441200 microseconds
        C solution from function avx256_row_matmul:
        Time taken by function: 225400 microseconds

    Where,

    constexpr int m = 70;
    constexpr int n = 90;
    constexpr int k = 40;

    std::unique_ptr OUTPUT:

        C matrix solution from function naive_row_matmul:
        Time taken by function: 1949797400 microseconds
        C solution from function avx256_row_matmul:
        Time taken by function: 878351600 microseconds

    Where,

    constexpr int m = 700;
    constexpr int n = 900;
    constexpr int k = 400;

    Solving matrix-matrix multiple: 
    C[m][k] = C[m][k] + A[m][n] * B[n][k]
*/

// #include <lapack.h>
// #include "lapack.h"
#include "matmul.hpp"
#include<chrono>
#include<memory>

template <typename T, std::size_t Rows, std::size_t Cols>
using Matrix = matrix::Matrix<T,Rows,Cols>;

constexpr int m = 70;
constexpr int n = 90;
constexpr int k = 40;

// Needed to find leading dimensional value at compile time
// Only useful for padding the arrays
// Would lead to some minor optimizations over handling fringe cases for tiled methods
#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

constexpr int mbuffer = (m % 4); 
constexpr int nbuffer = (n % 4); 
constexpr int kbuffer = (k % 4); 
constexpr int maxbuffer = max(mbuffer, max(nbuffer, kbuffer)); 
constexpr int ldmax = max(m + mbuffer, max(n + nbuffer, k + kbuffer));

int main() {
    //std::uniform_real_distribution<> dis = matrix::init_rand_gen(gen,0.0, 10.0);

    /* IGNORE BUFFERING MATRICES FOR NOW */
    // Matrix<double,ldmax,ldmax> C_buf = matrix::create_pad_matrix<double,m,k,ldmax>(C);
    // matrix::print_matrix<double,ldmax,ldmax>(C_buf);

    // std::cout << "ldmax = " << ldmax << std::endl;

    /*-------------------------------------------------------*/
    /*------------(double) std::array tests------------------*/
    /*-------------------------------------------------------*/

    Matrix<double,m,n> A; 
    matrix::fill_matrix<double,m,n>(1.0,A);
    // matrix::print_matrix<double,m,n>(A);
    // std::cout << std::endl;

    Matrix<double,n,k> B;
    matrix::fill_matrix<double,n,k>(1.0,B);
    // matrix::print_matrix<double,n,k>(B);
    // std::cout << std::endl;

    Matrix<double,m,k> C;
    matrix::fill_matrix<double,m,k>(0.0,C);
    // matrix::print_matrix<double,m,k>(C);
    // std::cout << std::endl;

    std::cout << "\nBeginning stack allocated std::array tests\n" << std::endl;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    naive_row_matmul<double>(m,n,k,A,B,C);
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

    /*-------------------------------------------------------*/
    /*-----------(double) array pointer tests----------------*/
    /*-------------------------------------------------------*/

    std::unique_ptr<Matrix<double,m,n>> A_ptr = 
    std::make_unique<Matrix<double,m,n>>();
    matrix::fill_matrix<double,m,n>(1.0,*(A_ptr));
    // matrix::print_matrix<double,m,n>(*(A_ptr));
    // std::cout << std::endl;

    std::unique_ptr<Matrix<double,n,k>> B_ptr = 
    std::make_unique<Matrix<double,n,k>>();
    matrix::fill_matrix<double,n,k>(1.0,*(B_ptr));
    // matrix::print_matrix<double,n,k>(*(B_ptr));
    // std::cout << std::endl;

    std::unique_ptr<Matrix<double,m,k>> C_ptr = 
    std::make_unique<Matrix<double,m,k>>();
    matrix::fill_matrix<double,m,k>(0.0,*(C_ptr));
    // matrix::print_matrix<double,m,k>(*(C_ptr));
    // std::cout << std::endl;

    std::cout << "\nBeginning heap std::array pointer tests\n" << std::endl;

    startTime = std::chrono::high_resolution_clock::now();
    naive_row_matmul<double>(m,n,k,*(A_ptr),*(B_ptr),*(C_ptr));
    endTime = std::chrono::high_resolution_clock::now();
    duration = endTime - startTime;

    std::cout << "C matrix solution from function naive_row_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
    // matrix::print_matrix<double,m,k>(*(C_ptr));

    matrix::fill_matrix<double,m,k>(0.0,*(C_ptr));

    startTime = std::chrono::high_resolution_clock::now();
    avx256_row_matmul(m,n,k,*(A_ptr),*(B_ptr),*(C_ptr));
    endTime = std::chrono::high_resolution_clock::now();
    duration = endTime - startTime;

    std::cout << "C solution from function avx256_row_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
    // matrix::print_matrix<double,m,k>(*(C_ptr));

    /*-------------------------------------------------------*/
    /*------------(float) std::array tests-------------------*/
    /*-------------------------------------------------------*/

    Matrix<float,m,n> A_flt; 
    matrix::fill_matrix<float,m,n>(1.0,A_flt);
    // matrix::print_matrix<double,m,n>(A);
    // std::cout << std::endl;

    Matrix<float,n,k> B_flt;
    matrix::fill_matrix<float,n,k>(1.0,B_flt);
    // matrix::print_matrix<double,n,k>(B);
    // std::cout << std::endl;

    Matrix<float,m,k> C_flt;
    matrix::fill_matrix<float,m,k>(0.0,C_flt);
    // matrix::print_matrix<double,m,k>(C);
    // std::cout << std::endl;

    std::cout << "\nBeginning <float> std::array tests\n" << std::endl;
    
    startTime = std::chrono::high_resolution_clock::now();
    naive_row_matmul<float>(m,n,k,A_flt,B_flt,C_flt);
    endTime = std::chrono::high_resolution_clock::now();
    duration = endTime - startTime;

    std::cout << "C matrix solution from function naive_row_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
    // matrix::print_matrix<double,m,k>(C);

    matrix::fill_matrix<float,m,k>(0.0,C_flt);

    startTime = std::chrono::high_resolution_clock::now();
    avx256_row_matmul(m,n,k,A_flt,B_flt,C_flt);
    endTime = std::chrono::high_resolution_clock::now();
    duration = endTime - startTime;

    std::cout << "C solution from function avx256_row_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
    // matrix::print_matrix<double,m,k>(C);
}
