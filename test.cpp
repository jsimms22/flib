// g++ test.cpp -lopenblas -lgfortran -mavx2

#include "matmul.hpp"
#include <openBLAS/cblas.h>
#include<chrono>
#include<memory>

template <typename T, std::size_t Rows, std::size_t Cols>
using Matrix = matrix::Matrix<T,Rows,Cols>;

constexpr int m = 4*256;//4;
constexpr int n = 8*256;//5;
constexpr int k = 3*256;//3;

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

    // Matrix<double,m,n> A; 
    // matrix::fill_matrix<double,m,n>(1.0,A);
    // // matrix::print_matrix<double,m,n>(A);
    // // std::cout << std::endl;

    // Matrix<double,n,k> B;
    // matrix::fill_matrix<double,n,k>(1.0,B);
    // // matrix::print_matrix<double,n,k>(B);
    // // std::cout << std::endl;

    // Matrix<double,m,k> C;
    // matrix::fill_matrix<double,m,k>(0.0,C);
    // // matrix::print_matrix<double,m,k>(C);
    // // std::cout << std::endl;

    // std::cout << "\nBeginning stack allocated std::array tests\n" << std::endl;
    
    // auto startTime = std::chrono::high_resolution_clock::now();
    // naive_row_matmul<double>(m,n,k,A,B,C);
    // auto endTime = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // std::cout << "C matrix solution from function naive_row_matmul:" << std::endl;
    // std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // // matrix::print_matrix<double,m,k>(C);

    // matrix::fill_matrix<double,m,k>(0.0,C);

    // startTime = std::chrono::high_resolution_clock::now();
    // avx256_row_matmul(m,n,k,A,B,C);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // std::cout << "C solution from function avx256_row_matmul:" << std::endl;
    // std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // // matrix::print_matrix<double,m,k>(C);

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

    auto startTime = std::chrono::high_resolution_clock::now();
    naive_row_matmul<double>(m,n,k,*(A_ptr),*(B_ptr),*(C_ptr));
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    std::cout << "C matrix solution from function naive_row_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // matrix::print_matrix<double,m,k>(*(C_ptr));

    matrix::fill_matrix<double,m,k>(0.0,*(C_ptr));

    startTime = std::chrono::high_resolution_clock::now();
    avx256_row_matmul(m,n,k,*(A_ptr),*(B_ptr),*(C_ptr));
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    std::cout << "C solution from function avx256_row_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // matrix::print_matrix<double,m,k>(*(C_ptr));

    /*-------------------------------------------------------*/
    /*------------(float) std::array tests-------------------*/
    /*-------------------------------------------------------*/

    // Matrix<float,m,n> A_flt; 
    // matrix::fill_matrix<float,m,n>(1.0,A_flt);
    // // matrix::print_matrix<float,m,n>(A_flt);
    // // std::cout << std::endl;

    // Matrix<float,n,k> B_flt;
    // matrix::fill_matrix<float,n,k>(1.0,B_flt);
    // // matrix::print_matrix<float,n,k>(B_flt);
    // // std::cout << std::endl;

    // Matrix<float,m,k> C_flt;
    // matrix::fill_matrix<float,m,k>(0.0,C_flt);
    // // matrix::print_matrix<float,m,k>(C_flt);
    // // std::cout << std::endl;

    // std::cout << "\nBeginning <float> std::array tests\n" << std::endl;
    
    // startTime = std::chrono::high_resolution_clock::now();
    // naive_row_matmul<float>(m,n,k,A_flt,B_flt,C_flt);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // std::cout << "C matrix solution from function naive_row_matmul:" << std::endl;
    // std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // // matrix::print_matrix<float,m,k>(C_flt);

    // matrix::fill_matrix<float,m,k>(0.0,C_flt);

    // startTime = std::chrono::high_resolution_clock::now();
    // avx256_row_matmul(m,n,k,A_flt,B_flt,C_flt);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // std::cout << "C solution from function avx256_row_matmul:" << std::endl;
    // std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // // matrix::print_matrix<float,m,k>(C_flt);

    /*-------------------------------------------------------*/
    /*-------------(int) std::array tests--------------------*/
    /*-------------------------------------------------------*/

    // Matrix<int,m,n> A_int; 
    // matrix::fill_matrix<int,m,n>(1.0,A_int);
    // // matrix::print_matrix<double,m,n>(A);
    // // std::cout << std::endl;

    // Matrix<int,n,k> B_int;
    // matrix::fill_matrix<int,n,k>(1.0,B_int);
    // // matrix::print_matrix<double,n,k>(B);
    // // std::cout << std::endl;

    // Matrix<int,m,k> C_int;
    // matrix::fill_matrix<int,m,k>(0.0,C_int);
    // // matrix::print_matrix<double,m,k>(C);
    // // std::cout << std::endl;

    // std::cout << "\nBeginning <int> std::array tests\n" << std::endl;
    
    // startTime = std::chrono::high_resolution_clock::now();
    // naive_row_matmul<int>(m,n,k,A_int,B_int,C_int);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // std::cout << "C matrix solution from function naive_row_matmul:" << std::endl;
    // std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // matrix::print_matrix<int,m,k>(C_int);

    // matrix::fill_matrix<int,m,k>(0,C_int);

    // startTime = std::chrono::high_resolution_clock::now();
    // avx256_row_matmul(m,n,k,A_int,B_int,C_int);
    // endTime = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    // std::cout << "C solution from function avx256_row_matmul:" << std::endl;
    // std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // matrix::print_matrix<int,m,k>(C_int);

    /*-------------------------------------------------------*/
    /*-----------------naive tile matmul---------------------*/
    /*-------------------------------------------------------*/

    // std::unique_ptr<Matrix<double,m,n>> A_ptr = 
    // std::make_unique<Matrix<double,m,n>>();
    // matrix::fill_matrix<double,m,n>(1.0,*(A_ptr));
    // // matrix::print_matrix<double,m,n>(*(A_ptr));
    // // std::cout << std::endl;

    // std::unique_ptr<Matrix<double,n,k>> B_ptr = 
    // std::make_unique<Matrix<double,n,k>>();
    // matrix::fill_matrix<double,n,k>(1.0,*(B_ptr));
    // // matrix::print_matrix<double,n,k>(*(B_ptr));
    // // std::cout << std::endl;

    std::unique_ptr<Matrix<double,m,k>> C_ptr_tile = 
    std::make_unique<Matrix<double,m,k>>();
    matrix::fill_matrix<double,m,k>(0.0,*(C_ptr_tile));
    // matrix::print_matrix<double,m,k>(*(C_ptr));
    // std::cout << std::endl;

    std::cout << "\nBeginning tile tests\n" << std::endl;
    
    startTime = std::chrono::high_resolution_clock::now();
    naive_tile_matmul<double>(m,n,k,*A_ptr,*B_ptr,*C_ptr_tile);
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    std::cout << "C matrix solution from function naive_tile_matmul:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;
    // matrix::print_matrix<float,m,k>(C_tile);
    // matrix::print_matrix<float,m,k>(C_flt);

    /*-------------------------------------------------------*/
    /*-------------------cblas dgemm-----------------------*/
    /*-------------------------------------------------------*/

    std::cout << "\nBeginning openBLAS dgemm reference benchmark\n" << std::endl;

    char TRANSA = 'N';
    char TRANSB = 'N';
    double ALPHA = 1.0;
    double BETA = 0;
    int LDA = m;
    int LDB = n;
    int LDC = k;

    double* A = new double[m*n];
    double* B = new double[n*k];
    double* C = new double[m*k];

    for (int i=0; i < (m*n); i++) { A[i] = i%3+1; }
    for (int i=0; i < (n*k); i++) { B[i] = i%3+1; }
    for (int i=0; i < (m*k); i++) { C[i] = 0; }

    startTime = std::chrono::high_resolution_clock::now();
    cblas_dgemm(CblasRowMajor,
                CblasTrans,
                CblasTrans, 
                m, k, k, 
                ALPHA, A, LDA, 
                B, LDB, BETA, 
                C, LDC);
    endTime = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    std::cout << "C matrix solution from function cblas_dgemm:" << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;

    delete[](A,B,C);
    return 0;
}
