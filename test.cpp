#include <iostream>
#include <cmath>
#include "dgemm.h"

using namespace std;

#define min(a,b) ((a < b))?(a):(b)
#define max(a,b) ((a > b))?(a):(b)

const int m = 2;//3
const int n = 3;//4
const int k = 1;//2

int main() {
  double** A;
  double** B;
  double** C;

  A = Build1(m,n,A);
  B = Build1(n,k,B);
  C = Build1(m,k,C); //this is simply to keep using Print3()

  Print3(m,n,k,A,B,C);

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

  cout << "mbuffer = " << mbuffer << "\n";
  cout << "nbuffer = " << nbuffer << "\n";
  cout << "kbuffer = " << kbuffer << "\n";
  cout << "maxbuffer = " << maxbuffer << "\n";
  cout << "ldmax = " << ldmax << "\n\n";

  cout << "Unrolled versions of the arrays\n\n";
  double* AA = Convert1(m,n,ldmax,A,true);
  cout << "A unrolled:\n";
  for(int i = 0; i < ldmax*ldmax; i++) {cout << AA[i] << " ";} cout << "\n\n";
  
  double* BB = Convert1(n,k,ldmax,B,true);
  cout<< "B unrolled:\n";
  for(int i = 0; i < ldmax*ldmax; i++) {cout << BB[i] << " ";} cout << "\n\n";

  double* CC = (double*)calloc((ldmax*ldmax),sizeof(double));
  cout << "C unrolled:\n";
  for(int i = 0; i < ldmax*ldmax; i++) {cout << CC[i] << " ";} cout << "\n\n";


  matmul(ldmax,m,n,k,AA,BB,CC);
  
  cout << "C solution from function matmul:\n";
  for(int i = 0; i < ldmax*ldmax; i++) {cout << CC[i] << " ";} cout << "\n\n";
  
  delete A,B,C,AA,BB,CC;
}
