#ifndef build3
#define build3

/* Header file for building 2d arrays, printing 2d arrays, 
 * and converting 2d to padded 1d arrays ---------------*/

#include<iostream>
#include<random>
#include<algorithm>
#include<array>

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(0.0, 10.0);

template <class T, std::size_t N>
void array_printer(int row, int column, std::array<T,N>& matrix) {
    for(int i = 0; i < column; i++) {
        for(int j = 0; j < row; j++) { 
            std::cout << matrix[i*column + j] << ' ';
        }
        std::cout << std::endl;
    }
}

template <class T, std::size_t N>
std::array<T,N> array_builder(double ALPHA) {
    std::array<T,N> matrix;
    for(int i = 0; i < matrix.size(); i++) { 
        matrix[i] = ALPHA * static_cast<T>(dis(gen)); 
    }
    return matrix;
}


double* Array_Buffer(int row, int column, int ldmax, double* array) {
  int i; int j;
  double* buffArray = (double*)calloc((ldmax*ldmax),sizeof(double));
    for(i = 0; i < row; i++) {
        for(j = 0; j < column; j++) {
            buffArray[i*ldmax + j] = array[i*column + j];
        }
    }
    return buffArray;
}

/*----- Function to find the leading dimensional value of 3 matrices -----*/
/*------ After buffering the 3 matrices to be equally sized squares ------*/
int ldmax_calc(int m, int n, int k) {
  int mbuffer = m % 4;
  mbuffer = (4 - mbuffer);
  int nbuffer = n % 4; 
  nbuffer = (4 - nbuffer);
  int kbuffer = k % 4;
  kbuffer = (4 - kbuffer);

  int maxbuffer = std::max(mbuffer,nbuffer);
  maxbuffer = std::max(maxbuffer,kbuffer);

  int ldmax = std::max(m+mbuffer,n+nbuffer); 
  ldmax = std::max(ldmax,k+kbuffer);
 
  printf("mbuffer = %d\n", mbuffer); 
  printf("nbuffer = %d\n", nbuffer);
  printf("kbuffer = %d\n", kbuffer);
  printf("maxbuffer = %d\n", maxbuffer);
  printf("ldmax = %d\n\n", ldmax);

  return ldmax;
}

#endif
