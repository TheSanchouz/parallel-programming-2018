#pragma once
#include <random>
#include <iostream>
#include <iomanip>
#include <omp.h>

class Matrix
{
private:
	int rows, cols;
	double **matrix;

	void Split4Strassen(Matrix &a11, Matrix &a12, Matrix &a21, Matrix &a22)
	{
		a11 = Submatrix(0, 0, rows / 2 - 1, cols / 2 - 1);
		a12 = Submatrix(0, cols / 2, rows / 2 - 1, cols - 1);
		a21 = Submatrix(rows / 2, 0, rows - 1, cols / 2 - 1);
		a22 = Submatrix(rows / 2, cols / 2, rows - 1, cols - 1);
	}
	static Matrix Collect4Strassen(const Matrix &a11, const Matrix &a12, const Matrix &a21, const Matrix &a22)
	{
		int rows = a11.rows;
		int cols = a11.cols;

		Matrix res(rows * 2, cols * 2);

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				res.matrix[i][j] = a11.matrix[i][j];
				res.matrix[i][j + cols] = a12.matrix[i][j];
				res.matrix[i + rows][j] = a21.matrix[i][j];
				res.matrix[i + rows][j + cols] = a22.matrix[i][j];
			}
		}

		return res;
	}

	void Peeling(Matrix &A11, Matrix &a12, Matrix &a21, Matrix &a22)
	{
		int oddRows, oddCols;

		if (rows % 2 == 0) oddRows = 0;
		else oddRows = 1;

		if (cols % 2 == 0) oddCols = 0;
		else oddCols = 1;

		A11 = Submatrix(0, 0, rows - 1 - oddRows, cols - 1 - oddCols);

		if (oddRows == 1 && oddCols == 1)
		{
			a12 = Submatrix(0, cols - 1, rows - 2, cols - 1);
			a21 = Submatrix(rows - 1, 0, rows - 1, cols - 2);
			a22 = Submatrix(rows - 1, cols - 1, rows - 1, cols - 1);
		}
		else if (oddRows == 1)
		{
			a12 = Matrix(rows - 1, 1);
			a21 = Submatrix(rows - 1, 0, rows - 1, cols - 1);
			a22 = Matrix(1, 1);
		}
		else if (oddCols == 1)
		{
			a12 = Submatrix(0, cols - 1, rows - 1, cols - 1);
			a21 = Matrix(1, cols - 1);
			a22 = Matrix(1, 1);
		}
		else
		{
			a12 = Matrix(rows, 1);
			a21 = Matrix(1, cols);
			a22 = Matrix(1, 1);
		}
	}
	static Matrix Fixup(Matrix &A11, const Matrix &a12, const Matrix &a21, const Matrix &a22,
		Matrix &B11, const Matrix &b12, const Matrix &b21, const Matrix &b22,
		const int rows, const int cols)
	{
		const Matrix C11 = _MultiplyStrassenEven(A11, B11) + a12 * b21;
		const Matrix c12 = A11 * b12 + a12 * b22;
		const Matrix c21 = a21 * B11 + a22 * b21;
		const Matrix c22 = a21 * b12 + a22 * b22;

		return CollectAfterFixup(C11, c12, c21, c22, rows, cols);
	}
	static Matrix CollectAfterFixup(const Matrix &C11, const Matrix &c12, const Matrix &c21, const Matrix &c22,
		const int rows, const int cols)
	{
		Matrix res(rows, cols);

		for (int i = 0; i < C11.rows; i++)
		{
			for (int j = 0; j < C11.cols; j++)
			{
				res.matrix[i][j] = C11.matrix[i][j];
			}
		}

		if (cols % 2 != 0)
		{
			for (int i = 0; i < c12.rows; i++)
			{
				res.matrix[i][res.cols - 1] = c12.matrix[i][0];
			}
		}

		if (rows % 2 != 0)
		{
			for (int i = 0; i < c21.cols; i++)
			{
				res.matrix[res.rows - 1][i] = c21.matrix[0][i];
			}
		}

		if (rows % 2 != 0 && cols % 2 != 0)
		{
			res.matrix[res.rows - 1][res.cols - 1] = c22.matrix[0][0];
		}

		return res;
	}

	static Matrix _MultiplyStrassen(Matrix &lhs, Matrix &rhs)
	{
		const static int THRESHOLD = 128 * 128;

		if (lhs.cols * rhs.rows <= THRESHOLD)
		{
			return lhs * rhs;
		}

		if (lhs.rows % 2 != 0 || lhs.cols % 2 != 0 || rhs.rows % 2 != 0 || rhs.cols % 2 != 0)
		{
			return _MultiplyStrassenOdd(lhs, rhs);
		}
		else
		{
			return _MultiplyStrassenEven(lhs, rhs);
		}
	}
	static Matrix _MultiplyStrassenEven(Matrix &lhs, Matrix &rhs)
	{
		int rowsA = lhs.rows / 2;
		int colsA = lhs.cols / 2;

		int rowsB = rhs.rows / 2;
		int colsB = rhs.cols / 2;

		Matrix a11(rowsA, colsA);
		Matrix a12(rowsA, colsA);
		Matrix a21(rowsA, colsA);
		Matrix a22(rowsA, colsA);

		Matrix b11(rowsB, colsB);
		Matrix b12(rowsB, colsB);
		Matrix b21(rowsB, colsB);
		Matrix b22(rowsB, colsB);

		lhs.Split4Strassen(a11, a12, a21, a22);
		rhs.Split4Strassen(b11, b12, b21, b22);

		Matrix t1 = a11 + a22;
		Matrix t2 = b11 + b22;
		Matrix t3 = a21 + a22;
		Matrix t4 = b12 - b22;
		Matrix t5 = b21 - b11;
		Matrix t6 = a11 + a12;
		Matrix t7 = a21 - a11;
		Matrix t8 = b11 + b12;
		Matrix t9 = a12 - a22;
		Matrix t10 = b21 + b22;

		Matrix p1, p2, p3, p4, p5, p6, p7;

		#pragma omp parallel sections
		{
			#pragma omp section
			{
				p1 = _MultiplyStrassen(t1, t2);
			}
			#pragma omp section
			{
				p2 = _MultiplyStrassen(t3, b11);
			}
			#pragma omp section
			{
				p3 = _MultiplyStrassen(a11, t4);
			}
			#pragma omp section
			{
				p4 = _MultiplyStrassen(a22, t5);
			}
			#pragma omp section
			{
				p5 = _MultiplyStrassen(t6, b22);
			}
			#pragma omp section
			{
				p6 = _MultiplyStrassen(t7, t8);
			}
			#pragma omp section
			{
				p7 = _MultiplyStrassen(t9, t10);
			}
		}

		Matrix c11 = p1 + p4 - p5 + p7;
		Matrix c12 = p3 + p5;
		Matrix c21 = p2 + p4;
		Matrix c22 = p1 - p2 + p3 + p6;

		return Collect4Strassen(c11, c12, c21, c22);
	}
	static Matrix _MultiplyStrassenOdd(Matrix &lhs, Matrix &rhs)
	{
		Matrix A11, a12, a21, a22;
		Matrix B11, b12, b21, b22;

		lhs.Peeling(A11, a12, a21, a22);
		rhs.Peeling(B11, b12, b21, b22);

		return Fixup(A11, a12, a21, a22, B11, b12, b21, b22, lhs.rows, rhs.cols);
	}

public:
	Matrix(int rows = 0, int cols = 0) : rows(rows), cols(cols)
	{
		matrix = new double*[rows];
		for (int i = 0; i < rows; i++)
		{
			matrix[i] = new double[cols] {0};
		}
	}
	Matrix(const Matrix &src)
	{
		rows = src.rows;
		cols = src.cols;
		matrix = new double*[rows];

		for (int i = 0; i < rows; i++)
		{
			matrix[i] = new double[cols];
			for (int j = 0; j < cols; j++)
			{
				matrix[i][j] = src.matrix[i][j];
			}
		}
	}
	Matrix &operator=(const Matrix &src)
	{
		if (rows != src.rows || cols != src.cols)
		{
			rows = src.rows;
			cols = src.cols;
			matrix = new double*[rows];
			for (int i = 0; i < rows; i++)
			{
				matrix[i] = new double[cols];
			}
		}

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				matrix[i][j] = src.matrix[i][j];
			}
		}

		return *this;
	}
	~Matrix()
	{
		for (int i = 0; i < rows; i++)
		{
			delete[] matrix[i];
		}
		delete[] matrix;
	}

	void InitRandom(int seed = 1)
	{
		static std::default_random_engine generator(seed);
		static std::uniform_int_distribution<> distribution(-10, 10);

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				matrix[i][j] = distribution(generator);
			}
		}
	}
	Matrix Transpose() const
	{
		Matrix res(cols, rows);

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				res.matrix[j][i] = matrix[i][j];
			}
		}

		return res;
	}
	Matrix Submatrix(int rowStart, int colStart, int rowEnd, int colEnd) const
	{
		Matrix res;

		if (rowEnd - rowStart >= 0 && colEnd - colStart >= 0)
		{
			res = Matrix(rowEnd - rowStart + 1, colEnd - colStart + 1);
			for (int i = rowStart; i <= rowEnd; i++)
			{
				for (int j = colStart; j <= colEnd; j++)
				{
					res.matrix[i - rowStart][j - colStart] = matrix[i][j];
				}
			}
		}

		return res;
	}

	friend bool operator==(const Matrix &lhs, const Matrix &rhs)
	{
		if (lhs.rows == rhs.rows && lhs.cols == rhs.cols)
		{
			for (int i = 0; i < lhs.rows; i++)
			{
				for (int j = 0; j < lhs.cols; j++)
				{
					if (lhs.matrix[i][j] != rhs.matrix[i][j])
					{
						return false;
					}
				}
			}

			return true;
		}

		return false;
	}

	friend Matrix operator+(const Matrix &lhs, const Matrix &rhs)
	{
		if (lhs.rows == rhs.rows && lhs.cols == rhs.cols)
		{
			Matrix res(lhs.rows, lhs.cols);

			for (int i = 0; i < res.rows; i++)
			{
				for (int j = 0; j < res.cols; j++)
				{
					res.matrix[i][j] = lhs.matrix[i][j] + rhs.matrix[i][j];
				}
			}

			return res;
		}
		else throw "Matrices are disproportionate";
	}
	friend Matrix operator-(const Matrix &lhs, const Matrix &rhs)
	{
		if (lhs.rows == rhs.rows && lhs.cols == rhs.cols)
		{
			Matrix res(lhs.rows, lhs.cols);

			for (int i = 0; i < res.rows; i++)
			{
				for (int j = 0; j < res.cols; j++)
				{
					res.matrix[i][j] = lhs.matrix[i][j] - rhs.matrix[i][j];
				}
			}

			return res;
		}
		else throw "Matrices are disproportionate";
	}
	friend Matrix operator*(const Matrix &lhs, const Matrix &rhs)
	{
		if (lhs.cols == rhs.rows)
		{
			Matrix res(lhs.rows, rhs.cols);

			#pragma omp parallel for
			for (int i = 0; i < res.rows; i++)
			{
				for (int j = 0; j < res.cols; j++)
				{
					for (int k = 0; k < lhs.cols; k++)
					{
						res.matrix[i][j] += lhs.matrix[i][k] * rhs.matrix[k][j];
					}
				}
			}

			return res;
		}
		else throw "Matrices are disproportionate";
	}
	friend Matrix MultiplyStrassen(Matrix &lhs, Matrix &rhs)
	{
		Matrix res;

		if (lhs.cols == rhs.rows)
		{
			res = _MultiplyStrassen(lhs, rhs);
		}

		return res;
	}

	friend std::ostream &operator<<(std::ostream &stream, const Matrix &src)
	{
		for (int i = 0; i < src.rows; i++)
		{
			for (int j = 0; j < src.cols; j++)
			{
				stream << std::setw(5) << src.matrix[i][j];
			}

			stream << std::endl;
		}

		return stream;
	}
};