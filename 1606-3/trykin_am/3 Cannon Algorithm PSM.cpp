// //3 Cannon Algorithm PSM.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.


#include "pch.h"
#include "mpi.h"
#include <random>
#include <iomanip>
#include <string>
#include <chrono>
#include <thread>

using namespace std;

inline double ScalarProduct(double *a, double *b, int size)
{
	double res = 0.0;

	for (int i = 0; i < size; i++)
	{
		res += a[i] * b[i];
	}

	return res;
}
void InitMatrix(double *matrix, int rows, int cols)
{
	std::default_random_engine generator(rand());
	std::uniform_int_distribution<> distribution(-10, 10);

	for (int i = 0; i < rows * cols; i++)
	{
		matrix[i] = distribution(generator);
	}
}
void TransposeMatrix(double *matrix, double *transMatrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			transMatrix[j * rows + i] = matrix[i * cols + j];
		}
	}
}
void PrintMatrix(double *matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::cout << setw(7) << matrix[i * cols + j] << " ";
		}
		std::cout << endl;
	}
}


int main(int argc, char *argv[])
{
	srand((unsigned int)time(nullptr));

	const int ROOT = 0;
	const int DATA_TAG_SEND_MATRIX_A_BLOCK = 1;
	const int DATA_TAG_SEND_MATRIX_B_BLOCK = 2;
	const int DATA_TAG_SEND_MATRIX_C_BLOCK = 3;
	const int DATA_TAG_SENDRECEIVE_MUTUAL_MATRIX_A_BLOCK = 4;
	const int DATA_TAG_SENDRECEIVE_MUTUAL_MATRIX_B_BLOCK = 5;

	int procNum;
	int procRank;
	double time;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	int orderMatrix = (argc == 1) ? 10 : stoi(argv[1]);

	double *A = nullptr;
	double *B = nullptr;
	double *C = nullptr;

	if (sqrt(procNum) != (int)sqrt(procNum))
	{
		if (procRank == ROOT)
		{
			cout << "Number of processes isn't a perfect squared!" << endl;
		}

		MPI_Finalize();
		return -1;
	}
	if (orderMatrix % (int)sqrt(procNum) != 0)
	{
		if (procRank == ROOT)
		{
			cout << "Blocks of matrices aren't squared!" << endl;
		}

		MPI_Finalize();
		return -2;
	}
	if (procRank == ROOT)
	{
		A = new double[orderMatrix * orderMatrix];
		B = new double[orderMatrix * orderMatrix];
		C = new double[orderMatrix * orderMatrix]{ 0 };

		InitMatrix(A, orderMatrix, orderMatrix);
		InitMatrix(B, orderMatrix, orderMatrix);

		bool showMatrix = false;
		if (argc == 3 && string(argv[2]) == "true")
		{
			showMatrix = true;
		}

		if (showMatrix)
		{
			cout << "Matrix A with size " << orderMatrix << "x" << orderMatrix << endl;
			PrintMatrix(A, orderMatrix, orderMatrix);
			cout << endl;

			cout << "Matrix B with size " << orderMatrix << "x" << orderMatrix << endl;
			PrintMatrix(B, orderMatrix, orderMatrix);
			cout << endl;
		}
	}

	MPI_Comm MPI_MY_COMM_CART;
	int dims[2] = { 0, 0 };
	const int periods[2] = { 1, 1 };

	MPI_Dims_create(procNum, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &MPI_MY_COMM_CART);

	int myCoords[2];
	MPI_Cart_coords(MPI_MY_COMM_CART, procRank, 2, myCoords);

	int blockSize = orderMatrix / (int)sqrt(procNum);

	double *bufA = new double[blockSize * blockSize]{ 0 };
	double *bufB = new double[blockSize * blockSize]{ 0 };
	double *bufC = new double[blockSize * blockSize]{ 0 };


	MPI_Datatype type;
	MPI_Datatype MPI_MY_DATATYPE_MATRIX_BLOCK;
	const int sizes[2]		= { orderMatrix, orderMatrix };
	const int subsizes[2]	= { blockSize, blockSize };
	const int starts[2]		= { 0, 0 };

	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
	MPI_Type_create_resized(type, 0, blockSize * sizeof(double), &MPI_MY_DATATYPE_MATRIX_BLOCK);
	MPI_Type_commit(&MPI_MY_DATATYPE_MATRIX_BLOCK);


	int *displasmentsA	= (procRank == ROOT) ? new int[procNum] : nullptr;
	int *displasmentsB	= (procRank == ROOT) ? new int[procNum] : nullptr;
	int *displasmentsC	= (procRank == ROOT) ? new int[procNum] : nullptr;
	int *sendCounts		= (procRank == ROOT) ? new int[procNum] : nullptr;

	time = MPI_Wtime();

	if (procRank == ROOT)
	{
		for (int i = ROOT; i < procNum; i++)
		{
			sendCounts[i] = 1;

			int coords[2];
			MPI_Cart_coords(MPI_MY_COMM_CART, i, 2, coords);

			int procForA;
			int procForB;
			MPI_Cart_rank(MPI_MY_COMM_CART, new int[2]{ coords[0], coords[1] + coords[0] }, &procForA);
			MPI_Cart_rank(MPI_MY_COMM_CART, new int[2]{ coords[0] + coords[1], coords[1] }, &procForB);

			int coordsA[2];
			int coordsB[2];
			MPI_Cart_coords(MPI_MY_COMM_CART, procForA, 2, coordsA);
			MPI_Cart_coords(MPI_MY_COMM_CART, procForB, 2, coordsB);

			displasmentsA[i] = coordsA[0] * orderMatrix + coordsA[1];
			displasmentsB[i] = coordsB[0] * orderMatrix + coordsB[1];
			displasmentsC[i] = coords[0] * orderMatrix + coords[1];
		}
	}

	MPI_Scatterv(A, sendCounts, displasmentsA, MPI_MY_DATATYPE_MATRIX_BLOCK, bufA, blockSize * blockSize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	MPI_Scatterv(B, sendCounts, displasmentsB, MPI_MY_DATATYPE_MATRIX_BLOCK, bufB, blockSize * blockSize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	//this_thread::sleep_for(chrono::milliseconds(500 * procRank));
	//std::cout << "process " << procRank << " with coord (" << myCoords[0] << ", " << myCoords[1] << ")"
	//	<< " has a submatrix A with size " << blockSize << "x" << blockSize << endl;
	//PrintMatrix(bufA, blockSize, blockSize);
	//std::cout << "process " << procRank << " with coord (" << myCoords[0] << ", " << myCoords[1] << ")"
	//	<< " has a submatrix B with size " << blockSize << "x" << blockSize << endl;
	//PrintMatrix(bufB, blockSize, blockSize);
	//std::cout << endl;

	for (int count = 0; count < (int)sqrt(procNum); count++)
	{
		double *transBufB = new double[blockSize * blockSize];
		TransposeMatrix(bufB, transBufB, blockSize, blockSize);

		for (int i = 0; i < blockSize; i++)
		{
			for (int j = 0; j < blockSize; j++)
			{
				bufC[i * blockSize + j] += ScalarProduct(bufA + i * blockSize, transBufB + j * blockSize, blockSize);
			}
		}

		//this_thread::sleep_for(chrono::milliseconds(400 * procRank));
		//std::cout << "process " << procRank << " with coord (" << myCoords[0] << ", " << myCoords[1] << ")"
		//	<< " has a submatrix A with size " << blockSize << "x" << blockSize << endl;
		//PrintMatrix(bufA, blockSize, blockSize);
		//std::cout << "process " << procRank << " with coord (" << myCoords[0] << ", " << myCoords[1] << ")"
		//	<< " has a submatrix B with size " << blockSize << "x" << blockSize << endl;
		//PrintMatrix(bufB, blockSize, blockSize);
		//std::cout << endl;


		int prevProc;
		int nextProc;

		MPI_Cart_shift(MPI_MY_COMM_CART, 1, -1, &prevProc, &nextProc);
		MPI_Status statusA;
		MPI_Sendrecv_replace(bufA, blockSize * blockSize, MPI_DOUBLE, nextProc, DATA_TAG_SENDRECEIVE_MUTUAL_MATRIX_A_BLOCK,
			prevProc, DATA_TAG_SENDRECEIVE_MUTUAL_MATRIX_A_BLOCK, MPI_COMM_WORLD, &statusA);

		MPI_Cart_shift(MPI_MY_COMM_CART, 0, -1, &prevProc, &nextProc);
		MPI_Status statusB;
		MPI_Sendrecv_replace(bufB, blockSize * blockSize, MPI_DOUBLE, nextProc, DATA_TAG_SENDRECEIVE_MUTUAL_MATRIX_B_BLOCK,
			prevProc, DATA_TAG_SENDRECEIVE_MUTUAL_MATRIX_B_BLOCK, MPI_COMM_WORLD, &statusB);
	}

	//this_thread::sleep_for(chrono::milliseconds(400 * procRank));
	//std::cout << "process " << procRank << " with coord (" << myCoords[0] << ", " << myCoords[1] << ")"
	//	<< " has a submatrix C with size " << blockSize << "x" << blockSize << endl;
	//PrintMatrix(bufC, blockSize, blockSize);
	//std::cout << endl;


	MPI_Gatherv(bufC, blockSize * blockSize, MPI_DOUBLE, C, sendCounts, displasmentsC, MPI_MY_DATATYPE_MATRIX_BLOCK, ROOT, MPI_COMM_WORLD);
	time = MPI_Wtime() - time;

	delete[] bufA;
	delete[] bufB;
	delete[] bufC;

	if (procRank == ROOT)
	{
		bool showMatrix = false;
		if (argc == 3 && string(argv[2]) == "true")
		{
			showMatrix = true;
		}

		if (showMatrix)
		{
			cout << "Matrix C with size " << orderMatrix << "x" << orderMatrix << endl;
			PrintMatrix(C, orderMatrix, orderMatrix);
		}

		std::cout.setf(std::ios::fixed);
		std::cout.precision(10);
		std::cout << "Time = " << time << std::endl;

		delete[] A;
		delete[] B;
		delete[] C;
	}

	MPI_Type_free(&MPI_MY_DATATYPE_MATRIX_BLOCK);
	MPI_Type_free(&type);

	MPI_Finalize();
	return 0;
}