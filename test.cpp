#include <iostream>
//#include "dgemm.h"
#include "build3.h"

using namespace std;

const int m = 3;
const int n = 4;
const int k = 2;

int main() {
  double** A;
  double** B;
  double** C;

  A = Build1(m,n,A);
  B = Build1(n,k,B);
  C = Build1(m,k,C);

  Print3(m,n,k,A,B,C);

  delete A,B,C;
}
