pragma once
#include <random>
#include <iomanip>
#include <iostream>

class Matrix
{
private:
	int rows, cols;
	double **matrix;

	void Split4Strassen(Matrix& a11, Matrix& a12, Matrix& a21, Matrix& a22)
	{
		a11 = Submatrix(0, 0, rows / 2 - 1, cols / 2 - 1);
		a12 = Submatrix(0, cols / 2, rows / 2 - 1, cols - 1);
		a21 = Submatrix(rows / 2, 0, rows - 1, cols / 2 - 1);
		a22 = Submatrix(rows / 2, cols / 2, rows - 1, cols - 1);
	}
	static Matrix Collect4Strassen(Matrix& a11, Matrix& a12, Matrix& a21, Matrix& a22)
	{
		int rows = a11.rows;
		int cols = a11.cols;

		Matrix res(rows << 1, cols << 1);

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
	static Matrix MultiplyStrassen(Matrix lhs, Matrix rhs, int size)
	{
		if (size <= 128)
		{
			return lhs * rhs;
		}

		size = size >> 1;

		Matrix a11(size, size);
		Matrix a12(size, size);
		Matrix a21(size, size);
		Matrix a22(size, size);

		Matrix b11(size, size);
		Matrix b12(size, size);
		Matrix b21(size, size);
		Matrix b22(size, size);

		lhs.Split4Strassen(a11, a12, a21, a22);
		rhs.Split4Strassen(b11, b12, b21, b22);

		Matrix p1 = MultiplyStrassen(a11 + a22, b11 + b22,	size);
		Matrix p2 = MultiplyStrassen(a21 + a22, b11,		size);
		Matrix p3 = MultiplyStrassen(a11,		b12 - b22,	size);
		Matrix p4 = MultiplyStrassen(a22,		b21 - b11,	size);
		Matrix p5 = MultiplyStrassen(a11 + a12, b22,		size);
		Matrix p6 = MultiplyStrassen(a21 - a11, b11 + b12,	size);
		Matrix p7 = MultiplyStrassen(a12 - a22, b21 + b22,	size);

		Matrix c11 = p1 + p4 - p5 + p7;
		Matrix c12 = p3 + p5;
		Matrix c21 = p2 + p4;
		Matrix c22 = p1 - p2 + p3 + p6;

		return Collect4Strassen(c11, c12, c21, c22);
	}
public:
	Matrix(int rows = 0, int cols = 0) : rows(rows), cols(cols)
	{
		matrix = new double*[rows];
		for (int i = 0; i < rows; i++)
		{
			matrix[i] = new double[cols];
		}
	}
	Matrix(const Matrix& src)
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
	Matrix& operator=(const Matrix& src)
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

	void InitRandom()
	{
		std::default_random_engine generator(rand());
		std::uniform_int_distribution<> distribution(-10, 10);

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				matrix[i][j] = distribution(generator);
			}
		}
	}

	double& operator()(int i, int j)
	{
		return matrix[i][j];
	}
	const double& operator()(int i, int j) const
	{
		return matrix[i][j];
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
	Matrix Submatrix(int rowStart, int colStart, int rowEnd, int colEnd)
	{
		Matrix res(rowEnd - rowStart + 1, colEnd - colStart + 1);

		for (int i = rowStart; i <= rowEnd; i++)
		{
			for (int j = colStart; j <= colEnd; j++)
			{
				res.matrix[i - rowStart][j - colStart] = matrix[i][j];
			}
		}

		return res;
	}

	Matrix operator+(const Matrix& rhs)
	{
		Matrix res;

		if (this->rows == rhs.rows && this->cols == rhs.cols)
		{
			res = Matrix(rows, cols);

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
					res.matrix[i][j] = this->matrix[i][j] + rhs.matrix[i][j];
				}
			}
		}

		return res;
	}
	Matrix operator-(const Matrix& rhs)
	{
		Matrix res;

		if (this->rows == rhs.rows && this->cols == rhs.cols)
		{
			res = Matrix(rows, cols);

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
					res.matrix[i][j] = this->matrix[i][j] - rhs.matrix[i][j];
				}
			}
		}

		return res;
	}
	Matrix operator*(const Matrix& rhs)
	{
		Matrix res;

		if (this->cols == rhs.rows)
		{
			res = Matrix(this->rows, rhs.cols);
			double *transposedCol = new double[rhs.rows];

			for (int j = 0; j < rhs.cols; j++)
			{
				for (int k = 0; k < this->cols; k++)
				{
					transposedCol[k] = rhs.matrix[k][j];
				}

				for (int i = 0; i < this->rows; i++)
				{
					double scalarProduct = 0;
					for (int k = 0; k < this->rows; k++)
					{
						scalarProduct += this->matrix[i][k] * transposedCol[k];
					}
					
					res.matrix[i][j] = scalarProduct;
				}
			}
		}

		return res;
	}
	friend Matrix MultiplyStrassen(Matrix& lhs, Matrix& rhs)
	{
		Matrix res;

		if (lhs.cols == rhs.rows)
		{
			res = MultiplyStrassen(lhs, rhs, lhs.cols);
		}

		return res;
	}

	friend std::ostream& operator<<(std::ostream& stream, const Matrix& src)
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