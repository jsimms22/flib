#include <stdio.h>
#include <float.h>
#include <math.h>
#include "cblas.h"
#include "matmul.h"
#include "build3.h"

#ifdef GETTIMEOFDAY
#include <sys/time.h>
#else
#include <time.h>
#endif


const int m = 3;
const int n = 4;
const int k = 2;
/*
#define cblas_dgemm dgemm
extern void dgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

void reference_dgemm(int N, double ALPHA, double* A, double* B, double* C) {
  enum CBLAS_ORDER ORDER = 102; //CBlasColMajor;
  enum CBLAS_TRANSPOSE TRANSA = 111; //CblasNoTrans;
  enum CBLAS_TRANSPOSE TRANSB = 111; //CblasNoTrans;
  int M = N;
  int K = N;
  double BETA = ALPHA;
  int LDA = N;
  int LDB = N;
  int LDC = N;
  dgemm (ORDER, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}*/

extern void matmul(int, int, int, int, double*, double*, double*);

int main() {
  double* A;
  double* B;
  double* C;
 
  A = Array_Builder(1.0,m,n);
  B = Array_Builder(1.0,n,k);
  
  /* below is used to validate the unrolled 2D array 
   * is printing to look like a 2D array correctly
  for(int i = 0; i < m*n; i++) {
    printf("%g ", A[i]);
  }
  printf("\n");
  printf("\n");
  */

  Array_Printer(m,n,A/*,1*/);
  
  /* below is used to validate the unrolled 2D array 
   * is printing to look like a 2D array correctly
  for(int i = 0; i < n*k; i++) {
    printf("%g ", B[i]);
  }
  printf("\n");
  printf("\n");
  */

  Array_Printer(n,k,B/*,1*/);

  int ldmax;
  ldmax = ldmax_calc(m,n,k);

  double* AA = Array_Buffer(m,n,ldmax,A/*,1*/);
  Array_Printer(ldmax,ldmax,AA/*,1*/);
  
  double* BB = Array_Buffer(n,k,ldmax,B/*,1*/);
  Array_Printer(ldmax,ldmax,BB/*,1*/);

  double* CC = (double*)calloc((ldmax*ldmax),sizeof(double));
  Array_Printer(ldmax,ldmax,CC/*,1*/);
  
  clock_t beginTime, endTime;

  beginTime = clock();
  matmul(ldmax,m,n,k,AA,BB,CC);
  endTime = clock();

  printf("matmul function time = %ld\n", endTime - beginTime);

  //CC = (double*)calloc((ldmax*ldmax),sizeof(double));

  /*beginTime = clock();
  reference_dgemm(n, -3.0*DBL_EPSILON*ldmax, AA, BB, CC);
  endTime = clock();
  */
  //printf("reference function time = %ld\n", endTime - beginTime);

  printf("C Solution from Function matmul: \n");
  Array_Printer(ldmax,ldmax,CC/*,1*/);

  //cout << "C solution from function matmul:\n";
  //for(int i = 0; i < ldmax*ldmax; i++) {cout << CC[i] << " ";} cout << "\n\n";
  
  free(A);
  free(B);
  free(C);
  free(AA);
  free(BB);
  free(CC);
}
