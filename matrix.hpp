#pragma once

/* Header file for building 2d arrays, printing 2d arrays, 
 * and converting 2d to padded 1d arrays ---------------*/

#include<iostream>
#include<random>
#include<algorithm>
#include<array>
#include<iomanip>

namespace matrix 
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 10.0); // Defining the desired random range

    // std::uniform_real_distribution<> init_rand_gen(std::mt19937 gen(), double lower, double upper) 
    // { 
    //     std::uniform_real_distribution<> dis(lower, upper); // Defining the desired random range
    //     return dis;
    // }

    template <class T, std::size_t Rows, std::size_t Cols>
    using Matrix = std::array<T, Rows*Cols>;

    template <class T, std::size_t Rows, std::size_t Cols>
    void print_matrix(Matrix<T,Rows,Cols> matrix) 
    {
        std::cout << "Printing unrolled matrix[" << Rows << "][" << Cols << "]" << std::endl;
        for(std::size_t i = 0; i < Rows; i++) {
            for(std::size_t j = 0; j < Cols; j++) { 
                std::cout << std::setw(10) << matrix[i*Cols + j] << ' ';
            }
            std::cout << std::endl;
        }
    }

    template <class T, std::size_t Rows, std::size_t Cols>
    void fill_matrix(double ALPHA, Matrix<T,Rows,Cols>& matrix) 
    {
        for(std::size_t i = 0; i < matrix.size(); i++) {
            // matrix[i] = i; // debugging print_matrix bounds
            matrix[i] = ALPHA * 1;// static_cast<T>(dis(gen)); 
        }
    }

    template <class T, std::size_t Rows, std::size_t Cols, std::size_t Ldmax>
    std::array<T,Ldmax*Ldmax> create_pad_matrix(Matrix<T,Rows,Cols>& matrix) 
    {
        // Create matrix of required buffered size
        std::array<T,Ldmax*Ldmax> buffered_matrix;
        // Initialize indices with 0 (0 will not affect matrix-matrix operations)
        fill_matrix<double,Ldmax,Ldmax>(0.0,buffered_matrix);

        for(std::size_t i = 0; i < Rows; i++) {
            for(std::size_t j = 0; j < Cols; j++) {
                buffered_matrix[i*Ldmax + j] = matrix[i*Cols + j];
            }
        }
        return buffered_matrix;
    }
}
