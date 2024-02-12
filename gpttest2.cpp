#include <iostream>
#include <immintrin.h>
#include <array>

// template <typename T, std::size_t Rows, std::size_t Cols>
// using Matrix = std::array<std::array<T, Cols>, Rows>;

template <typename T, std::size_t Rows, std::size_t Cols>
void matrixMultiplyAVX256(std::array<T, Cols*Rows>& A, std::array<T, Cols*Rows>& B, std::array<T, Cols*Rows>& C) {
    //constexpr std::size_t P = Rows; // P is the common dimension

    for (std::size_t i = 0; i < Rows*Cols; ++i) {
        for (std::size_t j = 0; j < Rows*Cols; j += 8) {
            __m256 result = _mm256_setzero_ps();

            for (std::size_t k = 0; k < Rows*Cols; ++k) {
                __m256 a = _mm256_broadcast_ss(&A[i+k*Rows]);
                __m256 b = _mm256_loadu_ps(&B[k+j*Rows]);
                result = _mm256_add_ps(result, _mm256_mul_ps(a, b));
            }

            _mm256_storeu_ps(&C[i+j*Rows], result);
        }
    }
}

int main() {
    // Example usage with float type
    constexpr std::size_t M = 3;
    constexpr std::size_t N = 4;
    //constexpr std::size_t P = 5;

    // Matrix<float, M, N> A;
    // Matrix<float, N, P> B;
    // Matrix<float, M, P> C;
    std::array<float, N*M> A;
    std::array<float, N*M> B;
    std::array<float, N*M> C;

    // Initialize matrices A and B with some values
    for (int i = 0; i < A.size(); i++) { A[i] = 1; }
    for (int i = 0; i < B.size(); i++) { B[i] = 1; }

    // Perform matrix multiplication using AVX-256
    matrixMultiplyAVX256<float,M,N>(A,B,C);

    // Display the result matrix C
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            std::cout << C[i+j*M] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
