#include <immintrin.h> // Header for AVX
#include <iostream>

// Function to perform matrix multiplication using AVX-256
void matrixMultiplyAVX256(float* A, float* B, float* C, int m, int n, int p) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; j += 8) { // Process 8 elements at a time with AVX-256
            __m256 result = _mm256_setzero_ps(); // Initialize result to zeros

            for (int k = 0; k < n; ++k) {
                // Load 8 elements from A[i][k] and broadcast to all elements in the AVX register
                __m256 a = _mm256_broadcast_ss(&A[i * n + k]);

                // Load 8 elements from B[k][j] and multiply with corresponding elements in 'a'
                __m256 b = _mm256_loadu_ps(&B[k * p + j]);

                // Multiply-add the elements
                result = _mm256_add_ps(result, _mm256_mul_ps(a, b));
            }

            // Store the result in the C matrix
            _mm256_storeu_ps(&C[i * p + j], result);
        }
    }
}

int main() {
    // Example usage
    const int m = 3; // Number of rows in matrix A
    const int n = 4; // Number of columns in matrix A and rows in matrix B
    const int p = 5; // Number of columns in matrix B

    // Allocate memory for matrices A, B, and C
    float A[m][n] = {{1.0, 2.0, 3.0, 4.0},
                     {5.0, 6.0, 7.0, 8.0},
                     {9.0, 10.0, 11.0, 12.0}};

    float B[n][p] = {{13.0, 14.0, 15.0, 16.0, 17.0},
                     {18.0, 19.0, 20.0, 21.0, 22.0},
                     {23.0, 24.0, 25.0, 26.0, 27.0},
                     {28.0, 29.0, 30.0, 31.0, 32.0}};

    float C[m][p];

    // Perform matrix multiplication using AVX-256
    matrixMultiplyAVX256(&A[0][0], &B[0][0], &C[0][0], m, n, p);

    // Display the result matrix C
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            std::cout << C[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
