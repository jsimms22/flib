#ifndef DGEMM
#define DGEMM

#include <vector>

#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

static inline void  2d_unroll(double** m, vector<double> v) {
  for(
}

static inline void naive_dgemm(int lda, int M, int N, int K, double* A, double* B, double* C) 
{ 
  for (int k = 0; k < lda; ++k) {
    for (int i = 0; i < lda; ++i) {
      register double cij = 0.0;
      for (int j = 0; j < lda; ++j) {
        cij += A[i] * B[j];
      }
      C[i] += cij;
    }
  }
}

// adding this to test lazygit
static inline void naive_dgemm(int lda, int M, int N, int K, double* A, double* B, double* C) { }

// Used for solving matmul C = C + A*B
// linear dimensional value will be n*m, where n and m are the 2-dim values

  /* for square: roll out, do method as built previously
   * for non-square: break into squares of largest possible size
   *                 roll out, do method as built previously
   * for fringe cases: we can either solve traditionally, or pad out fringe cases
   */

// lda = m*n, ldb = n*p, ldc = m*k
//
//
//    m*n        n*k        m*k
// |1 1 1 1|    |2 2|      |3 3|
// |1 1 1 1| *  |2 2|  =   |3 3|
// |1 1 1 1|    |2 2|      |3 3| 
//              |2 2|       
static inline void matmul(int m, int n, int k, double** A, double** B, double** C)
{
  register int BLOCK1 = 4; register int BLOCK2 = 8;
  register double zeroref = 0.0;


  int ldmax = max(m,n); ldmax = max(n,k); ldmax = max(m,k);
  
  int mbuffer = m%4; 
  int nbuffer = n%4; 
  int kbuffer = k%4;
    
  int maxbuffer = max(mbuffer,nbuffer); 
  maxbuffer = max(nbuffer,kbuffer);
  maxbuffer = max(mbuffer,kbuffer);

	/*-----------------------------------------------------------*/
	/*------------------NEED TO REWORK BUFFERING-----------------*/
	/*-----------------------------------------------------------*/
  
  vector<double> AA = malloc((ldmax + maxbuffer) * sizeof(zeroref));
  vector<double> BB = malloc((ldmax + maxbuffer) * sizeof(zeroref));
  vector<double> CC = malloc((ldmax + maxbuffer) * sizeof(zeroref));
  
  int counter = 0;
  int i = 0;
  for (; i < lda;) {
    AA[i] = A[i];
  }
  int i = 0;
  for (; i < ldb;) {
    BB[i] = B[i];
  }
  int i = 0;
  for (; i < ldc;) {
    CC[i] = C[i];
  }

  for (int x = 0; x < ldmax; x += BLOCK2) {
    int lim_i = x + min(BLOCK2,ldmax - x);
    for (int y = 0; y < ldmax; y += BLOCK2) {
      int lim_j = y + min(BLOCK2,ldmax - y);
      for (int z = 0; z < ldmax; z += BLOCK2) {
        int lim_k = z + min(BLOCK2,ldmax - z);
        for (int i = x; i < lim_i; i += BLOCK1) {
          int M = min(BLOCK1,lim_i - i);
          for (int j = y; j < lim_j; j += BLOCK1) {
            int N = min(BLOCK1,lim_j - j);
            for (int k = z; k < lim_k; k += BLOCK1) {
              int K = min(BLOCK1,lim_k - k);
              naive_dgemm(ldmax,K,N,M,&buffA[i + k*lda],&buffB[k + j*lda],&C[i + j*lda]);
            }
          }
        }
      }
    }
  }
}

#endif
