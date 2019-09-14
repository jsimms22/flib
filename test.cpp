#include <iostream>
#include "dgemm.h"

using namespace std;

main() {
	/*double matrix[5] = {0, 1, 2, 3, 4};
	double* A;
	A = &matrix[0];
	//double x = 10;

	//A = &x;

	cout << "*A = " << *A << "\n";
	cout << "A = " << A << "\n";

	cout << "*A + 1 = " << *(A+1) << "\n";
	cout << "A + 1 = " << (A + 1) << "\n";*/

	double A[6] = {1, 1, 1, 1, 1, 1};
	double B[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
	double C[6] = {0, 0, 0, 0, 0, 0};

	matmul(6,9,6,&A[0],&B[0],&C[0]);

	for (int i = 0; i < 16; ++i) {
		//for (int j = 1; j <= ; ++j) {
			cout << " " << C[i] << " ";
		//}	
		//cout << "\n";
	}

	return 0;
}
