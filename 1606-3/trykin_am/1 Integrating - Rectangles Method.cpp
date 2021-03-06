// 1 Integrating - Rectangles Method.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include "mpi.h"
#include <iostream>
#include <string>

double f(double x)
{
	return 4 / (1 + x * x);
}

double integral(double a, int n, double h)
{
	double integ = 0;

	for (int i = 0; i < n; i++)
	{
		integ += f(a + h * (0.5 + i)) * h;
	}

	return integ;
}

int main(int argc, char *argv[])
{
	int procNum, procRank;
	double time;

	double a = argc == 1 ? 0 : std::stod(argv[1]);
	double b = argc == 1 ? 1 : std::stod(argv[2]);
	int n = argc == 1 ? 1000000 : std::stoi(argv[3]);
	double result = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	time = MPI_Wtime();

	double h = (b - a) / n;
	double num = n / procNum;
	double localRange = (b - a) / procNum;
	double localA = a + procRank * localRange;
	double localResult = integral(localA, num, h);

	std::cout << "Process " << procRank << " has the partial result of " << 
		localResult << " and it's range in [" << localA << ", " << localA + localRange << "]" << std::endl;

	MPI_Reduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	time = MPI_Wtime() - time;

	if (procRank == 0)
	{
		std::cout.setf(std::ios::fixed);
		std::cout.precision(10);
		std::cout << "Integral = " << result << std::endl;
		std::cout << "Time - " << time;
	}

	MPI_Finalize();

	return 0;
}