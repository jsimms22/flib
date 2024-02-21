#pragma once

#include "matrix.hpp"
#include<mmintrin.h>  //MMX
#include<xmmintrin.h> //SSE
#include<emmintrin.h> //SSE2
#include<pmmintrin.h> //SSE3
#include<tmmintrin.h> //SSSE3
#include<smmintrin.h> //SSE4.1
#include<nmmintrin.h> //SSE4.2
#include<ammintrin.h> //SSE4A
#include<wmmintrin.h> //AES
#include<immintrin.h> //AVX, AVX2, FMA

 /*C Matrix 4x4     A Matrix   B Matrix
 * | 00 10 20 30 |  | 0x -> |  | 0x 1x 2x 3x |
 * | 01 11 21 31 |  | 1x -> |  |             |
 * | 02 12 22 32 |  | 2x -> |  |             |
 * | 03 13 23 33 |  | 3x -> |  |             |
 */
void do_4x4(int lda, int M, double* a, double* b, double* c) {
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
    
    for (int x = 0; x < M; ++x) {
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

void copy_a4 (int lda, const int M, std::array<double,16>& a_src, std::array<double,16>& a_dest) {
    for (int i = 0; i < M; i+=lda) {
        a_dest[i] = a_src[i];
        a_dest[i+1] = a_src[i+1];
        a_dest[i+2] = a_src[i+2];
        a_dest[i+3] = a_src[i+3];
        //a_src += lda;
    }
}

void copy_b4 (int lda, const int K, std::array<double,16>& b_src, std::array<double,16>& b_dest) {
    double *b_ptr0, *b_ptr1, *b_ptr2, *b_ptr3;
    b_ptr0 = b_src;
    b_ptr1 = b_ptr0 + lda;
    b_ptr2 = b_ptr1 + lda;
    b_ptr3 = b_ptr2 + lda;

    for (int i = 0; i < K; ++i) {
        b_dest++ = b_ptr0++;
        b_dest++ = b_ptr1++;
        b_dest++ = b_ptr2++;
        b_dest++ = b_ptr3++;
    }
}

// Double precision AVX256 row DGEMM
template <std::size_t MN, std::size_t NK, std::size_t MK>
void avx256_row_matmul(int M, int N, int K, std::array<double,MN>& A, std::array<double,NK>& B, std::array<double,MK>& C) 
{
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < K; j+=4) {
            __m256d result = _mm256_setzero_pd();  
            for (int k = 0; k < N; ++k) {
                __m256d a = _mm256_broadcast_sd(&A[k+i*N]);
                __m256d b = _mm256_loadu_pd(&B[j+k*K]);
                result = _mm256_add_pd(result, _mm256_mul_pd(a, b));
            }
            _mm256_storeu_pd(&C[j+i*K], result);
        }
    }
}

// floating point AVX256 row FGEMM
template <std::size_t MN, std::size_t NK, std::size_t MK>
void avx256_row_matmul(int M, int N, int K, std::array<float,MN>& A, std::array<float,NK>& B, std::array<float,MK>& C) 
{
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < K; j+=8) {
            __m256 result = _mm256_setzero_ps();  
            for (int k = 0; k < N; ++k) {
                __m256 a = _mm256_broadcast_ss(&A[k+i*N]);
                __m256 b = _mm256_loadu_ps(&B[j+k*K]);
                result = _mm256_add_ps(result, _mm256_mul_ps(a, b));
            }
            _mm256_storeu_ps(&C[j+i*K], result);
        }
    }
}

// Integer AVX256 row IGEMM
template <std::size_t MN, std::size_t NK, std::size_t MK>
void avx256_row_matmul(int M, int N, int K, std::array<int,MN>& A, std::array<int,NK>& B, std::array<int,MK>& C) 
{
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < K; ++j) {
            __m256i result = _mm256_setzero_si256();
            // __m256i result = _mm256_loadu_si256((__m256i*)&C[j+i*K]);
            for (int k = 0; k < N; k+=8) {
                __m256i a = _mm256_loadu_si256((__m256i*)&A[k+i*N]);
                __m256i b = _mm256_loadu_si256((__m256i*)&B[j+k*K]);
                result = _mm256_add_epi32(result, _mm256_mullo_epi32(a, b));
                // sum = _mm256_add_epi32(_mm256_loadu_si256((__m256i*)&C[j+i*K]), result);
            }
            _mm256_storeu_si256((__m256i*)&C[j+i*K], result);
            //C[j+i*K] = _mm256_add_epi32(_mm256_loadu_si256((__m256i*)&C[j+i*K]), result);
            //_mm256_storeu_si256((__m256i*)&C[j+i*K], _mm256_add_epi32(_mm256_loadu_si256((__m256i*)&C[j+i*K]), result));
        }
    }
}

// Used for solving matmul C = C + A*B
template <class T, std::size_t MN, std::size_t NK, std::size_t MK>
void naive_row_matmul(int M, int N, int K, std::array<T,MN>& A, std::array<T,NK>& B, std::array<T,MK>& C) 
{ 
    for (int i = 0; i < M; i++) {
        /* For each column j of B */
        for (int j = 0; j < K; j++) {
            /* Compute C(i,j) */
            double cji = 0;//C[j+i*N];
            for( int k = 0; k < N; k++) {
	            cji += A[k+i*N] * B[j+k*K];
            }
            C[j+i*K] = cji;
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
template <class T, std::size_t MN, std::size_t NK, std::size_t MK>
void matmul(int ldmax, int m, int n, int k, std::array<T,MN>& A, std::array<T,NK>& B, std::array<T,MK>& C)
{
    int BLOCK1 = 8; int BLOCK2 = 16;

    // row_naive_dgemm(ldmax,k,n,m,A,B,C);

    for (int x = 0; x < ldmax; x += BLOCK2) {
        int lim_i = x + std::min(BLOCK2,ldmax - x);
        for (int y = 0; y < ldmax; y += BLOCK2) {
            int lim_j = y + std::min(BLOCK2,ldmax - y);
            for (int z = 0; z < ldmax; z += BLOCK2) {
                int lim_k = z + std::min(BLOCK2,ldmax - z);
                for (int i = x; i < lim_i; i += BLOCK1) {
                    int M = std::min(BLOCK1,lim_i - i);
                    for (int j = y; j < lim_j; j += BLOCK1) {
                        int N = std::min(BLOCK1,lim_j - j);
                        for (int k = z; k < lim_k; k += BLOCK1) {
                            int K = std::min(BLOCK1,lim_k - k);
                            //blocked_row_dgemm(ldmax,M,N,K,A[k + i*ldmax],B[j + k*ldmax],C[j + i*ldmax]);
                        }
                    }
                }
            }
        }
    }
}