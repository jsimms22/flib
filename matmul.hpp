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

// Double precision AVX256 row DGEMM
// Used for solving matmul C = C + A*B
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

// Single precision AVX256 row FGEMM
// Used for solving matmul C = C + A*B
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
//
// template <class T, std::size_t MN, std::size_t NK, std::size_t MK>
// void matmul(int ldmax, int m, int n, int k, std::array<T,MN>& A, std::array<T,NK>& B, std::array<T,MK>& C)
// {
//     int BLOCK1 = 4; int BLOCK2 = 8;

//     row_naive_dgemm(ldmax,k,n,m,A,B,C);

//     for (int x = 0; x < ldmax; x += BLOCK2) {
//         int lim_i = x + std::min(BLOCK2,ldmax - x);
//         for (int y = 0; y < ldmax; y += BLOCK2) {
//             int lim_j = y + std::min(BLOCK2,ldmax - y);
//             for (int z = 0; z < ldmax; z += BLOCK2) {
//                 int lim_k = z + std::min(BLOCK2,ldmax - z);
//                 for (int i = x; i < lim_i; i += BLOCK1) {
//                     int M = std::min(BLOCK1,lim_i - i);
//                     for (int j = y; j < lim_j; j += BLOCK1) {
//                         int N = std::min(BLOCK1,lim_j - j);
//                         for (int k = z; k < lim_k; k += BLOCK1) {
//                             int K = std::min(BLOCK1,lim_k - k);
//                             //blocked_row_dgemm(ldmax,M,N,K,A[k + i*ldmax],B[j + k*ldmax],C[j + i*ldmax]);
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }