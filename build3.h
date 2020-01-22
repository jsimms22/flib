#ifndef build3
#define build3

/* Header file for building 2d arrays, printing 2d arrays, 
 * and converting 2d to padded 1d arrays ---------------*/

#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

/*--------- Simple function to print 3 2d arrays ------------*/
void Print3(int m, int n, int k, double** A, double**B, double**C) {
  int i; int j;

  cout << "\nA Matrix\n\n";
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      cout << A[i][j] << " ";
    }
    cout << "\n";
  }

  cout << "\nB Matrix\n\n";
  for(i = 0; i < n; i++) {
    for(j = 0; j < k; j++) {
      cout << B[i][j] << " ";
    }
    cout << "\n";
  }

  cout << "\nC Matrix\n\n";
  for(i = 0; i < m; i++) {
    for(j = 0; j < k; j++) {
      cout << C[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "\n";
}

/*------- Simple function to build and initialize 1 2d array --------*/
double** Build1(int row, int column, double** array) {
  int i; int j;

  srand(time(NULL));

  array = new double*[row];
  for(i = 0; i < row; i++) {
    array[i] = new double[column];
  }

  for(i = 0; i < row; i++) {
    for(j = 0; j < column; j++) {
      array[i][j] = (double)(rand() % 9 + 1);
    }
  }

  return array;
}

/*--------- Function to convert a 2d array into a 1d array ----------*/
//    *bool major determines if the 2d array will be 
//     unrolled column or row major
//    -> column = true/1
//    -> row = false/0
double* Convert1(int row, int column, int ldmax, double** array, bool major) {

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

#endif
