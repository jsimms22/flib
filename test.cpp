#include <iostream>
#include <stdio.h>
#include <string>

#include <float.h>
#include <math.h>
#include "build3.h"

//#include <lapack.h>
//#include "lapack.h"
#include "cblas.h"

#include "matmul.h"

#ifdef GETTIMEOFDAY
#include <sys/time.h>
#else
#include <time.h>
#endif

using namespace std;

#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

const int m = 3;//3
const int n = 4;//4
const int k = 2;//2

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

//extern const char* dgemm_desc;
extern void matmul(int, int, int, int, double*, double*, double*);

int main() {
  double** A;
  double** B;
  double** C;

  A = build1(1.0,m,n,A);
  B = build1(1.0,n,k,B);
  C = build1(1.0,m,k,C);

  print_2d_array(m,n,A);
  print_2d_array(n,k,B);
  print_2d_array(m,k,C);

  int ldmax = 0;
  ldmax_calc(m,n,k);

  //cout << "Unrolled versions of the arrays\n\n";
  double* AA = convert1(m,n,ldmax,A,true);
  //cout << "A unrolled:\n";
  //for(int i = 0; i < ldmax*ldmax; i++) {cout << AA[i] << " ";} cout << "\n\n";
  
  double* BB = convert1(n,k,ldmax,B,true);
  //cout<< "B unrolled:\n";
  //for(int i = 0; i < ldmax*ldmax; i++) {cout << BB[i] << " ";} cout << "\n\n";

  double* CC = (double*)calloc((ldmax*ldmax),sizeof(double));
  //cout << "C unrolled:\n";
  //for(int i = 0; i < ldmax*ldmax; i++) {cout << CC[i] << " ";} cout << "\n\n";

 // A = trans_mat(m,n,A);

  //print_2d_array(m,n,A);
  
  clock_t beginTime, endTime;
  
  matmul(ldmax,m,n,k,AA,BB,CC);

  beginTime = clock();
  matmul(ldmax,m,n,k,AA,BB,CC);
  endTime = clock();

  cout << "matmul function time = " << endTime - beginTime << "\n";

  //CC = (double*)calloc((ldmax*ldmax),sizeof(double));

  //beginTime = clock();
  //reference_dgemm(n, -3.0*DBL_EPSILON*ldmax, AA, BB, CC);
  //endTime = clock();
  
  //cout << "reference function time = " << endTime - beginTime << "\n";

  cout << "C solution from function matmul:\n";
  for(int i = 0; i < ldmax*ldmax; i++) {cout << CC[i] << " ";} cout << "\n\n";
  
  delete [] A,B,C,AA,BB,CC;
}
