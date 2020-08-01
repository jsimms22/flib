#ifndef build3
#define build3

/* Header file for building 2d arrays, printing 2d arrays, 
 * and converting 2d to padded 1d arrays ---------------*/

//#include <iostream>
//#include <cstdlib>
//#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

void Array_Printer(int row, int column, double* array, bool major) {
  int i; int j;

  for(i = 0; i < row; i++) {
    for(j = 0; j < column; j++) {
      if(major == 1) { printf("%g ", array[i*row + j]); }
      if(major == 0) { printf("%g ", array[j*row + i]); }
    }
    printf("\n");
  }
  printf("\n");
}

//--------- Simple function to print a 2d array ------------//
void print_2d_array(int row, int column, double** array) {
  int i, j;
  
  for(i = 0; i < row; i++) {
    for(j = 0; j < column; j++) {
      printf("%g ",array[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}


double* Array_Builder(double ALPHA, int row, int column) {
  int i;int j;
  double* array = (double*)malloc((row*column)*sizeof(double));

  srand(time(NULL));

  for(i = 0; i < row*column; i++) {
    array[i] = ALPHA * (double)(rand() % 9 + 1);
  }

  return array;
}

/*
//------- Simple function to build and initialize 1 2d array C++ Only --------//
double** build1(double ALPHA, int row, int column, double** array) {
  int i, j;

  srand(time(NULL));

  array = new double*[row];
  for(i = 0; i < row; i++) {
    array[i] = new double[column];
  }

  for(i = 0; i < row; i++) {
    for(j = 0; j < column; j++) {
      array[i][j] = ALPHA * (double)(rand() % 9 + 1);
    }
  }

  return array;
}*/


double* Array_Buffer(int row, int column, int ldmax, double* array, bool major) {
  int i; int j;
  double* buffArray = (double*)calloc((ldmax*ldmax),sizeof(double));

  if(major == 1) { //true means row major
    for(i = 0; i < ldmax; i++) {
      for(j = 0; j < ldmax; j++) {
        if(i < row && j < column) {
          buffArray[i*ldmax + j] = array[i*row + j];
        }
      }
    }
  }

  if(major == 0) { //false means column major
    for(i = 0; i < ldmax; i++) {
      for(j = 0; j < ldmax; j++) {
        if(j < row && i < column) {
          buffArray[ldmax*i + j] = array[j*row + i];
        }
      }
    }
  }

  return buffArray;
}

/*
//--------- Function to convert a 2d array into a 1d array C++ Only ----------//
//    *bool major determines if the 2d array will be 
//     unrolled column or row major
//    -> column = true/1
//    -> row = false
double* convert1(int row, int column, int ldmax, double** array, bool major) {

  int i; int j;

  double* unrolledArray = (double*)calloc((ldmax*ldmax),sizeof(double));

  if(major == false) {
    for(i = 0; i < ldmax; i++) {
      for(j = 0; j < ldmax; j++) {
        if(i < row && j < column) {
          unrolledArray[(ldmax)*i + j] = array[i][j];
        }
      }
    }
  }

  if(major == true) {
    for(i = 0; i < ldmax; i++) {
      for(j = 0; j < ldmax; j++) {
        if(j < row && i < column) {
          unrolledArray[(ldmax)*i + j] = array[j][i];
        }
      }
    }
  }

  return unrolledArray;
}
*/

/*----- Function to find the leading dimensional value of 3 matrices -----*/
/*------ After buffering the 3 matrices to be equally sized squares ------*/
int ldmax_calc(int m, int n, int k) {
  int mbuffer = m % 4;
  mbuffer = (4 - mbuffer);
  int nbuffer = n % 4; 
  nbuffer = (4 - nbuffer);
  int kbuffer = k % 4;
  kbuffer = (4 - kbuffer);

  int maxbuffer = max(mbuffer,nbuffer);
  maxbuffer = max(maxbuffer,kbuffer);

  int ldmax = max(m+mbuffer,n+nbuffer); 
  ldmax = max(ldmax,k+kbuffer);
 
  printf("mbuffer = %d\n", mbuffer); 
  printf("nbuffer = %d\n", nbuffer);
  printf("kbuffer = %d\n", kbuffer);
  printf("maxbuffer = %d\n", maxbuffer);
  printf("ldmax = %d\n\n", ldmax);

  return ldmax;
}


//results in a seg fault atm
/*double** trans_mat(int row, int column, double** array) {
  double** temp;
  build1(1.0,column,row,temp);

  int i, j;
  for(i = 0; i < row; i++) {
    for(j = 0; j < column; j++) {
      temp[j][i] = array[i][j];
    }
  }
  
  for(i = 0; i < column; i++) {
    for(j - 0; j < row; j++) {
      array[i][j] = temp[i][j];
    }
  }
  
  return temp;
}*/

#endif
