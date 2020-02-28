//#include <iostream>
#include <stdio.h>
//#include <string>

#include <float.h>
#include <math.h>

//#include <lapack.h>
//#include "lapacke.h"

#include "matmul.h"

#ifdef GETTIMEOFDAY
#include <sys/time.h>
#else
#include <time.h>
#endif

//using namespace std;

const int m = 3;//3
const int n = 4;//4
const int k = 2;//2

//#define DGEMM dgemm
/*extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

void reference_dgemm(int N, double ALPHA, double* A, double* B, double* C) {
  char TRANSA = 'N';
  char TRANSB = 'N';
  int M = N;
  int K = N;
  double BETA = 1.0;
  int LDA = N;
  int LDB = N;
  int LDC = N;
  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
}

extern const char* dgemm_desc;
extern void matmul(int, int, int, int, double*, double*, double*);*/

int main() {
  double* A;
  double* B;
  double* C;

  A = Array_Builder(1.0,m,n);
  B = Array_Builder(1.0,n,k);
  C = Array_Builder(1.0,m,k);

  Array_Printer(m,n,A,0);
  Array_Printer(n,k,B,0);
  Array_Printer(m,k,C,0);

  int ldmax;
  ldmax = ldmax_calc(m,n,k);

  //cout << "Unrolled versions of the arrays\n\n";
  double* AA = Array_Buffer(m,n,ldmax,A,0);
  //cout << "A unrolled:\n";
  //for(int i = 0; i < ldmax*ldmax; i++) {cout << AA[i] << " ";} cout << "\n\n";
  
  double* BB = Array_Buffer(n,k,ldmax,B,0);
  //cout<< "B unrolled:\n";
  //for(int i = 0; i < ldmax*ldmax; i++) {cout << BB[i] << " ";} cout << "\n\n";

  double* CC = (double*)calloc((ldmax*ldmax),sizeof(double));
  //cout << "C unrolled:\n";
  //for(int i = 0; i < ldmax*ldmax; i++) {cout << CC[i] << " ";} cout << "\n\n";

 // A = trans_mat(m,n,A);

  //print_2d_array(m,n,A);
  
  clock_t beginTime, endTime;
  
  //matmul(ldmax,m,n,k,AA,BB,CC);

  beginTime = clock();
  matmul(ldmax,m,n,k,AA,BB,CC);
  endTime = clock();

  printf("matmul function time = %ld\n", endTime - beginTime);

  /*CC = (double*)calloc((ldmax*ldmax),sizeof(double));

  beginTime = clock();
  reference_dgemm(n, -3.0*DBL_EPSILON*ldmax, AA, BB, CC);
  endTime = clock();
  
  printf("reference function time = %ld\n", endTime - beginTime);*/

  printf("C Solution from Function matmul: \n");
  Array_Printer(ldmax,ldmax,CC,0);

  //cout << "C solution from function matmul:\n";
  //for(int i = 0; i < ldmax*ldmax; i++) {cout << CC[i] << " ";} cout << "\n\n";
  
  //delete [] A,B,C,AA,BB,CC;
  free(A);
  free(B);
  free(C);
  free(AA);
  free(BB);
  free(CC);
}
