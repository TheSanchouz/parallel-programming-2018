#pragma once
#include "Matrix.h"
#include "tbb/task.h"

class StrassenMultiplier : public tbb::task
{
private:
	Matrix &A, &B, &C;
	int depth;

public:
	StrassenMultiplier(Matrix &A, Matrix &B, Matrix &C, int depth) :
		A(A), B(B), C(C), depth(depth) {}

	task* task::execute()
	{
		const static int THRESHOLD = 128 * 128;

		if (A.getCols() * B.getRows() <= THRESHOLD)
		{
			C = A * B;
		}
		else if (A.getRows() % 2 != 0 || A.getCols() % 2 != 0 || B.getRows() % 2 != 0 || B.getCols() % 2 != 0)
		{
			C = Matrix::_MultiplyStrassenOdd(A, B);
		}
		else
		{
			tbb::task_list list;
			int count = 1;

			int rowsA = A.getRows() / 2;
			int colsA = A.getCols() / 2;

			int rowsB = B.getRows() / 2;
			int colsB = B.getCols() / 2;

			Matrix a11(rowsA, colsA);
			Matrix a12(rowsA, colsA);
			Matrix a21(rowsA, colsA);
			Matrix a22(rowsA, colsA);

			Matrix b11(rowsB, colsB);
			Matrix b12(rowsB, colsB);
			Matrix b21(rowsB, colsB);
			Matrix b22(rowsB, colsB);

			A.Split4Strassen(a11, a12, a21, a22);
			B.Split4Strassen(b11, b12, b21, b22);

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

			// p1
			count++;
			list.push_back(*new (allocate_child())
				StrassenMultiplier(t1, t2, p1, depth + 1));

			// p2
			count++;
			list.push_back(*new (allocate_child())
				StrassenMultiplier(t3, b11, p2, depth + 1));

			// p3
			count++;
			list.push_back(*new (allocate_child())
				StrassenMultiplier(a11, t4, p3, depth + 1));

			// p4
			count++;
			list.push_back(*new (allocate_child())
				StrassenMultiplier(a22, t5, p4, depth + 1));

			// p5
			count++;
			list.push_back(*new (allocate_child())
				StrassenMultiplier(t6, b22, p5, depth + 1));

			// p6
			count++;
			list.push_back(*new (allocate_child())
				StrassenMultiplier(t7, t8, p6, depth + 1));

			// p7
			count++;
			list.push_back(*new (allocate_child())
				StrassenMultiplier(t9, t10, p7, depth + 1));

			set_ref_count(count);
			spawn_and_wait_for_all(list);

			Matrix c11 = p1 + p4 - p5 + p7;
			Matrix c12 = p3 + p5;
			Matrix c21 = p2 + p4;
			Matrix c22 = p1 - p2 + p3 + p6;

			C = Matrix::Collect4Strassen(c11, c12, c21, c22);
		}



		return NULL;
	}
};

void StrassenMultiplyTBB(Matrix &A, Matrix &B, Matrix &C)
{
	StrassenMultiplier& root = *new (tbb::task::allocate_root())
		StrassenMultiplier(A, B, C, 1);

	tbb::task::spawn_root_and_wait(root);
}