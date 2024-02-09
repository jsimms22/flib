#include "matmul.hpp"

// #ifdef GETTIMEOFDAY
// #include <sys/time.h>
// #else
// #include <time.h>
// #endif

// #define min(a,b) ((a < b))?(a):(b)
// #define max(a,b) ((a > b))?(a):(b)

constexpr int m = 2;//3
constexpr int n = 2;//4
constexpr int k = 2;//2

/*#define DGEMM dgemm_
extern "C" void dgemm_ (char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

void reference_dgemm(int N, double ALPHA, double* A, double* B, double* C) {
  char TRANSA = 'N';
  char TRANSB = 'N';
  int M = N;
  int K = N;
  double BETA = 1.0;
  int LDA = N;
  int LDB = N;
  int LDC = N;
  
  dgemm_ (&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
}*/

int main() {
    std::array<double,m*n> A = array_builder<double,m*n>(1.0);
    std::array<double,n*k> B = array_builder<double,n*k>(1.0);
    std::array<double,m*k> C;// = array_builder<double,m*k>(1.0);

    for (int i = 0; i < C.size(); i++) { C[i] = 0; }

    array_printer(m,n,A);
    array_printer(n,k,B);
    array_printer(m,k,C);

    int ldmax = ldmax_calc(m,n,k);

    //cout << "Unrolled versions of the arrays\n\n";
    //double* AA = convert1(m,n,ldmax,A,true);
    //cout << "A unrolled:\n";
    //for(int i = 0; i < ldmax*ldmax; i++) {cout << AA[i] << " ";} cout << "\n\n";
    
    //double* BB = convert1(n,k,ldmax,B,true);
    //cout<< "B unrolled:\n";
    //for(int i = 0; i < ldmax*ldmax; i++) {cout << BB[i] << " ";} cout << "\n\n";

    //double* CC = (double*)calloc((ldmax*ldmax),sizeof(double));
    //cout << "C unrolled:\n";
    //for(int i = 0; i < ldmax*ldmax; i++) {cout << CC[i] << " ";} cout << "\n\n";

    // A = trans_mat(m,n,A);

    //print_2d_array(m,n,A);
    
    //clock_t beginTime, endTime;
    
    matmul<double,m*n,n*k,m*k>(2,m,n,k,A,B,C);

    // beginTime = clock();
    // matmul(ldmax,m,n,k,A,B,C);
    // endTime = clock();

    // std::cout << "matmul function time = " << endTime - beginTime << "\n";

    //CC = (double*)calloc((ldmax*ldmax),sizeof(double));

    //beginTime = clock();
    //reference_dgemm(n, -3.0*DBL_EPSILON*ldmax, AA, BB, CC);
    //endTime = clock();
    
    //cout << "reference function time = " << endTime - beginTime << "\n";

    std::cout << "C solution from function matmul:\n";
    for (int i = 0; i < m*k; i++) {
        std::cout << C[i] << " ";
    } 
    std::cout << "\n\n";
}
