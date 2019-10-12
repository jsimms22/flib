#ifndef build3
#define build3

#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

void Print3(int m, int n, int k, double** A, double**B, double**C) {
  int i; int j;

  cout << "\nA Matrix\n\n";
  for(i = 0; i < m; i++) {
    cout << A[i][0] << " " << A[i][1] << " " << A[i][2] << " " << A[i][3] << "\n";
  }

  cout << "\nB Matrix\n\n";
  for(i = 0; i < n; i++) {
    cout << B[i][0] << " " << B[i][1] << "\n";
  }

  cout << "\nC Matrix\n\n";
  for(i = 0; i < m; i++) {
    cout << C[i][0] << " " << C[i][1] << "\n";
  }
  cout << "\n";
}

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

#endif
