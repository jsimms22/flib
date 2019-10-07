#ifndef DGEMM
#define DGEMM

#include <vector>

#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

// adding this to test lazygit
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

// Used for solving matmul C = C + A*B
// linear dimensional value will be n*m, where n and m are the 2-dim values

static inline void matmul(int lda, int ldb, int ldc, double* A, double* B, double* C)
{
  int BLOCK1 = 4;//256;
  int BLOCK2 = 8;//512;
  double zeroref = 0.0;

  int ldmax = max(lda,ldb);
  ldmax = max(ldb,ldc);
  ldmax = max(lda,ldc);
  
  int abuffer = lda%4;
  int bbuffer = ldb%4;
  int cbuffer = ldc%4;
    
  int maxbuffer = max(abuffer,bbuffer);
  maxbuffer = max(bbuffer,cbuffer);
  maxbuffer = max(abuffer,cbuffer);
  
  double buffA[ldmax + maxbuffer];
  double buffB[ldmax + maxbuffer];
  double buffC[ldmax + maxbuffer];

	printf("test 1\n");

	/*-----------------------------------------------------------*/
	/*------------------NEED TO REWORK BUFFERING-----------------*/
	/*-----------------------------------------------------------*/

  vector<double> AA = calloc((ldmax + maxbuffer), sizeof(zeroref));
  vector<double> BB = calloc((ldmax + maxbuffer), sizeof(zeroref));
  vector<double> CC = calloc((ldmax + maxbuffer), sizeof(zeroref));
  
  int counter = 0;
  int i = 0;
  for (; i < ldmax+maxbuffer;) {
    AA[i] = A
  }
  
  /*------------------old method for padding below-------------*/

  if(lda != ldb && lda != ldc) {
    int counter = 0;
    int i = 0;
    if (abuffer != 0 || lda < ldmax) {
      for(; counter < ldmax;)
      {
	buffA[i] = *(A + counter);
	i++; counter++;
	if (ldmax%i != 0) {
	  for (int j = 0; j < abuffer; j++) { 	
	    buffA[i] = zeroref; // this may be slow, not contiguous memory
	    i++;
	  }
	}
      }
    }
    if (bbuffer != 0 || ldb < ldmax) {
      i = 0;
      for(; counter < ldmax;)
      {
	buffB[i] = *(B + counter);
	i++; counter++;
	if (ldmax%i != 0) {
	  for (int j = 0; j < bbuffer; j++) { 	
	    buffB[i] = zeroref; // this may be slow, not contiguous memory
	    i++;
	  }
	}
      }
    }

    if (cbuffer != 0 || ldc < ldmax) {
      i = 0;
      for(; counter < ldmax;)
      {
	buffB[i] = *(B + counter);
	i++; counter++;
	if (ldmax%i != 0) {
	  for (int j = 0; j < bbuffer; j++) { 	
	    buffB[i] = zeroref; // this may be slow, not contiguous memory
	    i++;
	  }
	}
      }
    }
  }

	printf("test 2\n");

  /*for (int x = 0; x < ldmax; x += BLOCK2) {
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
  }*/
}

#endif
