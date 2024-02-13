#include <iostream>
// #include "matrix.hpp"
// #include <lapack.h>
// #include "lapack.h"
#include "matmul.hpp"

// #ifdef GETTIMEOFDAY
// #include <sys/time.h>
// #else
// #include <time.h>
// #endif

template <typename T, std::size_t Rows, std::size_t Cols>
using Matrix = matrix::Matrix<T,Rows,Cols>;

constexpr int m = 4;//3
constexpr int n = 4;//4
constexpr int k = 3;//2

// Needed to find leading dimensional value at compile time
// Only useful for padding the arrays
// Would lead to some minor optimizations over handling fringe cases
#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

constexpr int mbuffer = 4 - (m % 4); 
constexpr int nbuffer = 4 - (n % 4); 
constexpr int kbuffer = 4 - (k % 4); 
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

    matrix::print_matrix<double,m,n>(A);
    std::cout << std::endl;
    matrix::print_matrix<double,n,k>(B);
    std::cout << std::endl;
    matrix::print_matrix<double,m,k>(C);
    std::cout << std::endl;

    /* IGNORE BUFFERING MATRICES FOR NOW */
    // Matrix<double,ldmax,ldmax> C_buf = matrix::create_pad_matrix<double,m,k,ldmax>(C);
    // matrix::print_matrix<double,ldmax,ldmax>(C_buf);
    
    row_naive_dgemm<double,m*n,n*k,m*k>(ldmax,m,n,k,A,B,C);
    std::cout << "C solution from function matmul:" << std::endl;
    matrix::print_matrix<double,m,k>(C);
}
