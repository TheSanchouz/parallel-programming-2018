#include "pch.h"
#include "Matrix.h"
#include <iostream>

using namespace std;

int main()
{
	int l, m, n;
	int num_threads = 4;
	
	cout << "Enter a sizes of matrixes" << endl;
	cin >> l >> m >> n;

	cout << "You will multiply matrices with sizes " << l << "x" << m << " and " << m << "x" << n << endl;

	cout << "Set a number of threads" << endl;
	cin >> num_threads;
	omp_set_num_threads(num_threads);

	Matrix A(l, m);
	Matrix B(m, n);
	A.InitRandom();
	B.InitRandom();

	//Matrix C = A * B;
	double time = omp_get_wtime();
	Matrix D = MultiplyStrassen(A, B);
	cout << omp_get_wtime() - time << endl;

	//if (C == D)
	{
		cout << "\x1B[32m" << "OK! Multiply Strassen is correct!" << "\033[0m\t\t" << endl;
	}
	//else
	{
		cout << "\x1B[31m" << "FAILED! Multiply Strassen is uncorrect!" << "\033[0m\t\t" << endl;
	}

	return 0;
}