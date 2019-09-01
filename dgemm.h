#ifndef DGEMM
#define DGEMM

#define min(a,b) ((a < b))?(a):(b)

// adding this to test lazygit command
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

// Used for solving matmul D = C + A*B
// linear dimensional value will be n*m, where n and m are the 2-dim values

static inline void matmul(int lda, int ldb, int ldc, double* A, double* B, double* C)
{
  register int BLOCK1 = 256;
  register int BLOCK2 = 512;

  if(lda != ldb && lda != ldc) {
    int buffa = lda%2;
    int buffb = ldb%2;
    int buffc = ldc%2;

    if (buffa != 0) {}

    if (buffb != 0) {}

    if (buffa != 0) {}
  }

  for (int x = 0; x < lda; x + BLOCK2) {
    for (int y = 0; y < lda; y + BLOCK2) {
      for (int z = 0; z < lda; z + BLOCK2) {
        int lim_i = min(BLOCK2,lda - x);
        int lim_j = min(BLOCK2,lda - y);
        int lim_k = min(BLOCK2,lda - z);
        for (int i = 0; i < lim_i; i + BLOCK1) {
          for (int j = 0; j < lim_j; j + BLOCK1) {
            for (int k = 0; k < lim_k; k + BLOCK1) {
              int M = min(BLOCK1,lim_i - i);
              int N = min(BLOCK1,lim_j - j);
              int K = min(BLOCK1,lim_k - k);

              naive_dgemm(lda,M,N,K,A + i + k*lda,B + k + j*lda,C + i + j*lda);
            }
          }
        }
      }
    }
  }
}

#endif
