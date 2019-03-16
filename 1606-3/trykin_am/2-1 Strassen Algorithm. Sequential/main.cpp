#include "pch.h"
#include <omp.h>
#include "Matrix.h"
#include <iostream>
#include <iomanip>

using namespace std;


int main()
{
	srand((unsigned int)time(nullptr));

	int size;
	cin >> size;

	Matrix A(size, size);
	Matrix B(size, size);
	Matrix C(size, size);
	A.InitRandom();
	B.InitRandom();

	//cout << A << endl;
	//cout << B << endl;

	double time = omp_get_wtime();
	C = MultiplyStrassen(A, B);
	time = omp_get_wtime() - time;

	//cout << C << endl;
	cout << time << endl;

	return 0;
}