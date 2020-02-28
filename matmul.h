#ifndef DGEMM
#define DGEMM

#include "build3.h"
#include <mmintrin.h>
#include <immintrin.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <emmintrin.h>
#include <string.h>

const char* dgemm_desc = "AVX and Blocked DGEMM Function.";

#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

#define turn_even(x) (((x) & 1) ? (x+1) : (x))

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#define ARRAY(A,i,j) (A)[(j)*lda + (i)]

 /*C Matrix 4x4     A Matrix   B Matrix
 * | 00 10 20 30 |  | 0x -> |  | 0x 1x 2x 3x |
 * | 01 11 21 31 |  | 1x -> |  |             |
 * | 02 12 22 32 |  | 2x -> |  |             |
 * | 03 13 23 33 |  | 3x -> |  |             |
 */
static void do_4x4 (int lda, int K, double* a, double* b, double* c) {
  register __m256d a0x_3x,
    bx0, bx1, bx2, bx3,
    c00_30, c01_31,
    c02_32, c03_33;
  
  double* c01_31_ptr = c + lda;
  double* c02_32_ptr = c01_31_ptr + lda;
  double* c03_33_ptr = c02_32_ptr + lda;
  
  c00_30 = _mm256_loadu_pd(c);
  c01_31 = _mm256_loadu_pd(c01_31_ptr);
  c02_32 = _mm256_loadu_pd(c02_32_ptr);
  c03_33 = _mm256_loadu_pd(c03_33_ptr);
  
  for (int x = 0; x < K; ++x) {
    a0x_3x = _mm256_loadu_pd(a);
    a += 4;

    bx0 = _mm256_broadcast_sd(b++);
    bx1 = _mm256_broadcast_sd(b++);
    bx2 = _mm256_broadcast_sd(b++);
    bx3 = _mm256_broadcast_sd(b++);

    c00_30 = _mm256_add_pd(c00_30, _mm256_mul_pd(a0x_3x,bx0));
    c01_31 = _mm256_add_pd(c01_31, _mm256_mul_pd(a0x_3x,bx1));
    c02_32 = _mm256_add_pd(c02_32, _mm256_mul_pd(a0x_3x,bx2));
    c03_33 = _mm256_add_pd(c03_33, _mm256_mul_pd(a0x_3x,bx3));
  }
  
  _mm256_storeu_pd(c,c00_30);
  _mm256_storeu_pd(c01_31_ptr,c01_31);
  _mm256_storeu_pd(c02_32_ptr,c02_32);
  _mm256_storeu_pd(c03_33_ptr,c03_33);
}

static inline void copy_a4 (int lda, const int K, double* a_src, double* a_dest) {
  for (int i = 0; i < K; ++i) {
    *a_dest++ = *a_src;
    *a_dest++ = *(a_src + 1);
    *a_dest++ = *(a_src + 2);
    *a_dest++ = *(a_src + 3);
    a_src += lda;
  }
}

static inline void copy_b4 (int lda, const int K, double* b_src, double* b_dest) {
  double *b_ptr0, *b_ptr1, *b_ptr2, *b_ptr3;
  b_ptr0 = b_src;
  b_ptr1 = b_ptr0 + lda;
  b_ptr2 = b_ptr1 + lda;
  b_ptr3 = b_ptr2 + lda;

  for (int i = 0; i < K; ++i) {
    *b_dest++ = *b_ptr0++;
    *b_dest++ = *b_ptr1++;
    *b_dest++ = *b_ptr2++;
    *b_dest++ = *b_ptr3++;
  }
}

static void do_avx256 (int lda, int M, int N, int K, double* a, double* b, double* c) {
  __m256d m0,m1,m2,m3;
  for (int i = 0; i < M; i += 4) {
    for (int j = 0; j < N; ++j) {
      m0 = _mm256_setzero_pd();  
      for (int k = 0; k < K; ++k) {
	m1 = _mm256_load_pd(a+i+k*lda);
	m2 = _mm256_broadcast_sd(b+k+j*lda);
	m3 = _mm256_mul_pd(m1,m2);
	m0 = _mm256_add_pd(m0,m3);
      }
      m1 = _mm256_load_pd(c+i+j*lda);
      m0 = _mm256_add_pd(m0,m1);
      _mm256_storeu_pd(c+i+j*lda,m0);
    }
  }
}

static void inline blocked_column_dgemm(int lda, int M, int N, int K, double* A, double* B, double* C)
{
  double A_block[M*K], B_block[K*N];
  double *a_ptr, *b_ptr, *c;

/*
 * 8x8 blocks: slower than 4x4 blocks
 */
  /*int Nmax = N-7;
  int Mmax = M-7;
  int fringe1 = M%8;
  int fringe2 = N%8;
  
  int i = 0, j = 0, p = 0;
    
  for (j = 0; j < Nmax; j += 8) {
    b_ptr = &B_block[j*K];
    copy_b8 (lda, K, B + j*lda, b_ptr);
    for (i = 0; i < Mmax; i += 8) {
      a_ptr = &A_block[i*K];
      if (j == 0) copy_a8 (lda, K, A + i, a_ptr);
      c = C + i + j*lda;
      do_8x8 (lda, K, a_ptr, b_ptr, c);
    }
  }*/

 /* 4x4 blocks */
  int Nmax = N-3;
  int Mmax = M-3;
  int fringe1 = M%4;
  int fringe2 = N%4;

  register int i = 0, j = 0, p = 0;

  for (j = 0 ; j < Nmax; j += 4) 
  {
    b_ptr = &B_block[j*K];
    copy_b4(lda, K, B + j*lda, b_ptr);
    for (i = 0; i < Mmax; i += 4) {
      a_ptr = &A_block[i*K];
      if (j == 0) copy_a4(lda, K, A + i, a_ptr);
      c = C + i + j*lda;
      do_4x4(lda, K, a_ptr, b_ptr, c);
    }
  }

 /* 2x2 blocks */
  /*int Nmax = N-1;
  int Mmax = M-1;
  int fringe1 = M%2;
  int fringe2 = N%2;
  
  int i = 0, j = 0, p = 0;
  
  for (j = 0; j < Nmax; j += 2) {
    b_ptr = &B_block[j*K];
    copy_b2 (lda, K, B + j*lda, b_ptr);
    for (i = 0; i < Mmax; i += 2) {
      a_ptr = &A_block[i*K];
      if (j == 0) copy_a2 (lda, K, A + i, a_ptr);
      c = C + i + j*lda;
      do_2x2 (lda, K, a_ptr, b_ptr, c);
    }
  }*/

  /* Handle "fringes" */
  if (fringe1 != 0) {
    /* For each row of A */
    for ( ; i < M; ++i)
      /* For each column of B */ 
      for (p = 0; p < N; ++p) {
        /* Compute C[i,j] */
        register double c_ip = ARRAY(C,i,p);
        for (int k = 0; k < K; ++k)
          c_ip += ARRAY(A,i,k) * ARRAY(B,k,p);
        ARRAY(C,i,p) = c_ip;
      }
  }

  if (fringe2 != 0) {
    Mmax = M - fringe1;
    /* For each column of B */
    for ( ; j < N; ++j)
      /* For each row of A */ 
      for (i = 0; i < Mmax; ++i) {
        /* Compute C[i,j] */
        register double cij = ARRAY(C,i,j);
        for (int k = 0; k < K; ++k)
          cij += ARRAY(A,i,k) * ARRAY(B,k,j);
        ARRAY(C,i,j) = cij;
      }
  }
}

static inline void blocked_column_naive_dgemm(int lda, int M, int N, int K, double* A, double* B, double* C) {
  // For each row of A
  for (int i = 0; i < M; ++i) {
    // For each column of B
    for (int j = 0; j < N; ++j) {
      // Compute C[i,j] 
      register double cij = 0.0;
      for (int k = 0; k < K; ++k){
        cij += A[i+k*lda] * B[k+j*lda];
      }
      C[i+j*lda] += cij;
    }
  }
}

static inline void column_naive_dgemm(int lda, int M, int N, int K, double* A, double* B, double* C) 
{ 
  for (int i = 0; i < lda; i++) {
    /* For each column j of B */
    for (int j = 0; j < lda; j++) {
      /* Compute C(i,j) */
      double cij = C[i+j*lda];
        for( int k = 0; k < lda; k++) {
	        cij += A[i+k*lda] * B[k+j*lda];
        }
      C[i+j*lda] = cij;
    }
  }
}

//Below solves for row major DGEMM
static inline void row_naive_dgemm(int lda, int M, int N, int K, double* A, double* B, double* C) 
{ 
  for (int i = 0; i < lda; i++) {
    /* For each column j of B */
    for (int j = 0; j < lda; j++) {
      /* Compute C(i,j) */
      double cji = C[j+i*lda];
        for( int k = 0; k < lda; k++) {
	        cji += A[k+i*lda] * B[j+k*lda];
        }
      C[j+i*lda] = cji;
    }
  }
}

// Used for solving matmul C = C + A*B
// linear dimensional value will be n*m, where n and m are the 2-dim values

  /* for square: roll out, do method as built previously
   * for non-square: break into squares of largest possible size
   *                 roll out, do method as built previously
   * for fringe cases: we can either solve traditionally, or pad out fringe cases
   */

// lda = m*n, ldb = n*p, ldc = m*k
//
//    m*n        n*k        m*k
// |1 1 1 1|    |2 2|      |3 3|
// |1 1 1 1| *  |2 2|  =   |3 3|
// |1 1 1 1|    |2 2|      |3 3| 
//              |2 2|       
//
static inline void matmul(int ldmax, int m, int n, int k, double* AA, double* BB, double* CC)
{
  register int BLOCK1 = 4; register int BLOCK2 = 8;
  
  /*int mbuffer = m % 4;
  mbuffer = (4 - mbuffer);
  int nbuffer = n % 4;
  nbuffer = (4 - nbuffer);
  int kbuffer = k % 4;
  kbuffer = (4 - kbuffer);

  int maxbuffer = max(mbuffer,nbuffer);
  maxbuffer = max(maxbuffer,kbuffer);

  int ldmax = max(m+mbuffer,n+nbuffer); 
  ldmax = max(ldmax,k+kbuffer);*/
  
  //row_naive_dgemm(ldmax,k,n,m,AA,BB,CC);
  //column_naive_dgemm(ldmax,k,n,m,AA,BB,CC);

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
              blocked_column_dgemm(ldmax,M,N,K,&AA[i + k*ldmax],&BB[k + j*ldmax],&CC[i + j*ldmax]);
            }
          }
        }
      }
    }
  }
}

#endif
