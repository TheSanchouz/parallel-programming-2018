#include "pch.h"
#include <iostream>
#include <random>
#include <iomanip>
#include <string>
#include "mpi.h"

int mod(int a, int b)
{
	if (b < 0) return mod(a, -b);
	int ret = a % b;
	if (ret < 0)
		ret += b;
	return ret;
}

int scalarProduct(int *a, int *b, int size)
{
	int res = 0;

	for (int i = 0; i < size; i++)
	{
		res += a[i] * b[i];
	}

	return res;
}

void initMatrix(int *matrix, int rows, int cols)
{
	std::default_random_engine generator(rand());
	std::uniform_int_distribution<> distribution(-10, 10);

	for (int i = 0; i < rows * cols; i++)
	{
		matrix[i] = distribution(generator);
	}
}

void transposeMatrix(int *matrix, int *transMatrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			transMatrix[j * rows + i] = matrix[i * cols + j];
		}
	}
}

void transposeMatrix(int *matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::swap(matrix[i * cols + j], matrix[j * rows + i]);
		}
	}
}

void printMatrix(int *matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::cout << std::setw(5) << matrix[i * cols + j] << " ";
		}

		std::cout << std::endl;
	}
}

int main(int argc, char *argv[])
{
	srand((unsigned int)time(nullptr));

	const int root = 0;
	const int DATA_TAG = 0;
	double time;

	MPI_Status status;
	int procRank;
	int procNum;

	int rowsA = argc == 1 ? 3 : std::stoi(argv[1]);
	int colsA_rowsB = argc == 1 ? 4 : std::stoi(argv[2]);
	int colsB = argc == 1 ? 5 : std::stoi(argv[3]);

	int *matrixA = nullptr;
	int *matrixB = nullptr;
	int *transMatrixB = nullptr;
	int *matrixC = nullptr;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);


	if (procNum > rowsA)
	{
		if (procRank == root)
		{
			std::cout << "Error! Enter the number of process less than the number of rows of matrix A!" << std::endl;
		}

		MPI_Finalize();
		return -4;
	}

	if (procNum > colsB)
	{
		if (procRank == root)
		{
			std::cout << "Error! Enter the number of process less than the number of cols of matrix B!" << std::endl;
		}

		MPI_Finalize();
		return -4;
	}

	if (procRank == root)
	{
		//выделение памяти под матрицы
		/*matrixA = new int[rowsA * colsA_rowsB]{ 5, 7, 7, 2,
												7, 0, 1, 1,
												0, 5, 4, 6 };
		matrixB = new int[colsA_rowsB * colsB]{ 6, 10,  6,  6,  4,
												2,  0,  3,  5,  4,
												9,  1,  2,  8,  8,
												10, 2,  9,  4,  5 };
		transMatrixB = new int[colsB * colsA_rowsB];
		matrixC = new int[rowsA * colsB];*/

		matrixA = new int[rowsA * colsA_rowsB];
		matrixB = new int[colsA_rowsB * colsB];
		transMatrixB = new int[colsB * colsA_rowsB];
		matrixC = new int[rowsA * colsB];

		//заполнение матриц случайными числами
		initMatrix(matrixA, rowsA, colsA_rowsB);
		initMatrix(matrixB, colsA_rowsB, colsB);

		transposeMatrix(matrixB, transMatrixB, colsA_rowsB, colsB);

		bool showMatrix = false;
		if (argc == 5 && std::string(argv[4]) == "true")
		{
			showMatrix = true;
		}

		//вывод матриц на экран
		if (showMatrix)
		{
			std::cout << "Matrix A with size " << rowsA << "x" << colsA_rowsB << std::endl;
			printMatrix(matrixA, rowsA, colsA_rowsB);
			std::cout << std::endl;

			std::cout << "Matrix B with size " << colsA_rowsB << "x" << colsB << std::endl;
			printMatrix(matrixB, colsA_rowsB, colsB);
			std::cout << std::endl;

			//printMatrix(transMatrixB, colsB, colsA_rowsB);
		}
		else
		{
			//std::cout << "Matrix A with size " << rowsA << "x" << colsA_rowsB << std::endl;
			//std::cout << "Matrix B with size " << colsA_rowsB << "x" << colsB << std::endl;
		}
	}

	time = MPI_Wtime();

	//определение кол-ва строк и столбцов, по крайней мере минимальное, каждого процесса
	int everyHasRows = rowsA / procNum;
	int everyHasCols = colsB / procNum;

	//определение кол-ва строк и столбцов, которые необходимо дополнительно распределить к процессам
	int additiveRows = rowsA % procNum;
	int additiveCols = colsB % procNum;

	//массивы, содержащие кол-во посылаемых элементов для каждого процесса
	int *sizeOfSendRowsElem = procRank == root ? new int[procNum] : nullptr;
	int *sizeOfSendColsElem = procRank == root ? new int[procNum] : nullptr;

	//массивы, содержащие смещение
	int *displasmentsRowsElem = procRank == root ? new int[procNum] : nullptr;
	int *displasmentsColsElem = procRank == root ? new int[procNum] : nullptr;

	int *sizeOfReceiveRowsElem = procRank == root ? new int[procNum] : nullptr;
	int *displasmentsReceive = procRank == root ? new int[procNum] : nullptr;

	//мастер-процесс определяет эти массивы
	if (procRank == root)
	{
		//std::cout << "EveryHasRows " << everyHasRows << std::endl;
		//std::cout << "EveryHasCols " << everyHasCols << std::endl;
		//std::cout << "AdditiveRows " << additiveRows << std::endl;
		//std::cout << "AdditiveCols " << additiveCols << std::endl;
		//std::cout << std::endl;

		for (int i = 0; i < procNum; i++)
		{
			//определение числа передаваемых элеметов для процесса
			//если ранг процесса меньше чем число доп.строк и/или столбцов, то добавим для текущего процесса одну доп.строку и/или столбец
			sizeOfSendRowsElem[i] = i < additiveRows ? everyHasRows + 1 : everyHasRows;
			sizeOfSendColsElem[i] = i < additiveCols ? everyHasCols + 1 : everyHasCols;
			//определение числа присылаемых элементов на главный процесс
			sizeOfReceiveRowsElem[i] = sizeOfSendRowsElem[i] * colsB;

			//умножаем кол-во строк и столбцов на размер одной строки и одного столбца соотвественно
			//таким образом определяя кол-во передаваемых элементов процессу
			sizeOfSendRowsElem[i] *= colsA_rowsB;
			sizeOfSendColsElem[i] *= colsA_rowsB;

			//определение смещений
			displasmentsRowsElem[i] = i == 0 ? 0 : displasmentsRowsElem[i - 1] + sizeOfSendRowsElem[i - 1];
			displasmentsColsElem[i] = i == 0 ? 0 : displasmentsColsElem[i - 1] + sizeOfSendColsElem[i - 1];
			//определение смещения для итогового сбора
			displasmentsReceive[i] = i == 0 ? 0 : displasmentsReceive[i - 1] + sizeOfReceiveRowsElem[i - 1];

			//std::cout << "For process - " << i << ":" << std::endl;
			//std::cout << "	countOfRows			" << sizeOfSendRowsElem[i] / colsA_rowsB << std::endl;
			//std::cout << "	countOfCols			" << sizeOfSendColsElem[i] / colsA_rowsB << std::endl;
			//std::cout << "	sizeOfRowsElem			" << sizeOfSendRowsElem[i] << std::endl;
			//std::cout << "	sizeOfColsElem			" << sizeOfSendColsElem[i] << std::endl;
			//std::cout << "	displasmentsRow			" << displasmentsRowsElem[i] << std::endl;
			//std::cout << "	displasmentsCol			" << displasmentsColsElem[i] << std::endl;
			//std::cout << "	countOfReceiveRowsElem		" << sizeOfReceiveRowsElem[i] << std::endl;
			//std::cout << "	displasmentsReceive		" << displasmentsReceive[i] << std::endl;
			//std::cout << std::endl;
		}
	}

	//time = MPI_Wtime();

	//определение числа получаемых строк и столбцов для процесса
	//если ранг процесса меньше чем число доп.строк и/или столбцов, то добавим для текущего процесса одну доп.строку и/или столбец
	int receiveCountRows = procRank < additiveRows ? everyHasRows + 1 : everyHasRows;
	int receiveCountCols = procRank < additiveCols ? everyHasCols + 1 : everyHasCols;

	int *bufA = new int[receiveCountRows * colsA_rowsB]{ 0 };
	int *bufB = new int[receiveCountCols * colsA_rowsB]{ 0 };
	int *bufC = new int[receiveCountRows * colsB]{ 0 };

	//std::cout << "I am process - " << procRank << ", bufC[" << receiveCountRows * colsB << "]" << std::endl;

	//рассылка строк и столбцов
	MPI_Scatterv(matrixA, sizeOfSendRowsElem, displasmentsRowsElem, MPI_INT, bufA, receiveCountRows * colsA_rowsB, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Scatterv(transMatrixB, sizeOfSendColsElem, displasmentsColsElem, MPI_INT, bufB, receiveCountCols * colsA_rowsB, MPI_INT, root, MPI_COMM_WORLD);

	//определение предыдущего и следующего процесса
	int prevProc = mod(procRank - 1, procNum);
	int nextProc = mod(procRank + 1, procNum);

	/*std::cout << "I am process - " << procRank << ", for me" << std::endl
		<< "prevProc - " << prevProc << ", nextProc - " << nextProc << std::endl;*/

	double timeCalc = MPI_Wtime();
	for (int count = 0; count < procNum; count++)
	{
		int offsetCols = 0;
		for (int i = 0; i < mod(procRank - count, procNum); i++)
		{
			offsetCols += i < additiveCols ? everyHasCols + 1 : everyHasCols;
		}

		//std::cout << "For process " << procRank << ", count = " << count << ", offsetCols = " << offsetCols << std::endl;

		//скалярное произведение
		for (int i = 0; i < receiveCountRows; i++)
		{
			for (int j = 0; j < receiveCountCols; j++)
			{
				bufC[i * colsB + j + offsetCols] = scalarProduct(bufA + i * colsA_rowsB, bufB + j * colsA_rowsB, colsA_rowsB);
				//std::cout << "For process " << procRank << " [" << i * colsB << ", " << j + offsetCols << "]" 
				//	<< " result = " << scalarProduct(bufA + i * colsA_rowsB, bufB + j * colsA_rowsB, colsA_rowsB) << std::endl;
			}
		}
		//номера предыдущего и следующего процессов используются для совмещенного приема и передачи сообщения соответственно
		//топология - кольцо, передача и прием других столбцов

		//обмен должен происходить только (procNum - 1) раз, на последней итерации его не должно быть
		//обмен происходит между вычислением скалярного произведения
		if (count < procNum - 1)
		{
 			int newReceiveCountCols = mod(prevProc - count, procNum) < additiveCols ? everyHasCols + 1 : everyHasCols;

			//std::cout << "for process " << procRank << ", and count = " << count << ", newRCC = " << newReceiveCountCols << std::endl;

			//выделение нового буфера для приема сообщения
			int *newReceiveBuf = new int[newReceiveCountCols * colsA_rowsB];

			MPI_Sendrecv(bufB, receiveCountCols * colsA_rowsB, MPI_INT, nextProc, DATA_TAG,
				newReceiveBuf, newReceiveCountCols * colsA_rowsB, MPI_INT, prevProc, DATA_TAG, 
				MPI_COMM_WORLD, &status);

			delete bufB;
			bufB = new int[newReceiveCountCols * colsA_rowsB];
			memcpy(bufB, newReceiveBuf, newReceiveCountCols * colsA_rowsB * sizeof(int));

			receiveCountCols = newReceiveCountCols;
		}
	}

	//сбор строк результирующей матрицы
	MPI_Gatherv(bufC, receiveCountRows * colsB, MPI_INT, matrixC, sizeOfReceiveRowsElem, displasmentsReceive, MPI_INT, root, MPI_COMM_WORLD);

	time = MPI_Wtime() - time;

	delete[] bufA;
	delete[] bufB;
	delete[] bufC;
	std::cout << "Time calc = " << timeCalc << std::endl;
	if (procRank == root)
	{
		bool showMatrix = false;
		if (argc == 5 && std::string(argv[4]) == "true")
		{
			showMatrix = true;
		}

		std::cout.setf(std::ios::fixed);
		std::cout.precision(10);
		std::cout << "Time = " << time << std::endl;
		//std::cout << "Matrix C = A * B with size " << rowsA << "x" << colsB << std::endl;

		if (showMatrix)
		{
			printMatrix(matrixC, rowsA, colsB);
		}

		delete[] matrixA;
		delete[] matrixB;
		delete[] transMatrixB;
		delete[] matrixC;
	}

	MPI_Finalize();

	return 0;
}