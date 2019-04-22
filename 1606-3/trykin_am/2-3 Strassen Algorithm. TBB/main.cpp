#include "pch.h"
#include "Matrix.h"
#include "StrassenMultiplier.h"
#include <iostream>

#include "tbb/task_scheduler_init.h"
#include "tbb/tick_count.h"

using namespace std;

int main()
{
	int l, m, n;
	cout << "Enter a sizes of matrixes" << endl;
	cin >> l >> m >> n;

	cout << "You will multiply matrices with sizes " << l << "x" << m << " and " << m << "x" << n << endl;

	int num_threads;
	cout << "Set a number of threads" << endl;
	cin >> num_threads;

	Matrix A(l, m);
	Matrix B(m, n);
	A.InitRandom();
	B.InitRandom();

	//Matrix C = A * B;
	Matrix D;

	tbb::task_scheduler_init init(num_threads);
	tbb::tick_count parallel1 = tbb::tick_count::now();
	StrassenMultiplyTBB(A, B, D);
	tbb::tick_count parallel2 = tbb::tick_count::now();
	cout << "Time parallel = " << (parallel2 - parallel1).seconds() << " sec" << endl;

	//if (C == D)
	//{
	//	cout << "\x1B[32m" << "OK! Multiply Strassen is correct!" << "\033[0m\t\t" << endl;
	//}
	//else
	//{
	//	cout << "\x1B[31m" << "FAILED! Multiply Strassen is uncorrect!" << "\033[0m\t\t" << endl;
	//}

	return 0;
}