#include "matrix_inversion_parallel.h"
#include "helpers/common.h"
#include <omp.h>

// #include "/opt/homebrew/opt/libomp/include/omp.h" for MAC as this import is not working as usual

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* getcwd */
#include <stdbool.h>
#include <dirent.h>	  /* readdir */
#include <string.h>	  /* strlen */
#include <math.h>	  /* fabs */
#include <sys/time.h> /* gettimeofday */

/* Function for benchmarking */
void benchmark_matrix_inversion(int nrow, int ncol, double mat[nrow][ncol])
{
	struct timeval start, end;
	double mat_inv[nrow][ncol];

	/* Ensure threads are synchronized before timing */
#pragma omp parallel
	{
#pragma omp barrier
	}

	/* Start timing */
	gettimeofday(&start, NULL);

	/* Perform matrix inversion */
	if (!invert_matrix_par(nrow, ncol, mat, mat_inv))
	{
		printf("Matrix inversion failed during benchmarking.\n");
		return;
	}

	/* Ensure threads are synchronized after execution */
#pragma omp parallel
	{
#pragma omp barrier
	}

	/* End timing */
	gettimeofday(&end, NULL);

	/* Calculate elapsed time in milliseconds */
	double elapsed_time = (end.tv_sec - start.tv_sec) * 1000.0; /* Seconds to milliseconds */
	elapsed_time += (end.tv_usec - start.tv_usec) / 1000.0;		/* Microseconds to milliseconds */

	printf("Matrix inversion (Parallel) completed in %.3f ms for %dx%d matrix.\n", elapsed_time, nrow, ncol);
}

bool invert_matrix_par(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol])
{
	int n = nrow;

	/* Augment identity */
	double mat_aug[n][2 * ncol];
	augment_mat_par(n, mat, mat_aug);

	/* Forward elimination */
	bool res = gaussian_elimination_par(n, 2 * n, mat_aug);
	if (!res)
	{
		printf("GE failed\n");
		return false;
	}

	bool res2 = rref_par(n, 2 * n, mat_aug);
	if (!res2)
	{
		printf("RREF failed");
		return false;
	}

	extract_inverse_par(n, 2 * n, mat_aug, mat_inv);

	return true;
}

/* Extract the inverse part of the augmented matrix: parallelized by rows */
void extract_inverse_par(int nrow, int ncol, double mat_aug[nrow][ncol], double mat_inv[nrow][nrow])
{
	int i, j;
	for (i = 0; i < nrow; i++)
	{
#pragma omp parallel for
		for (j = 0; j < ncol; j++)
		{
			/* Copy only the right side of the augmented matrix */
			mat_inv[i][j] = mat_aug[i][nrow + j];
		}
	}
}

/* Normalize pivots and clear nonzero values below diagonal */
bool gaussian_elimination_par(int nrow, int ncol, double mat[nrow][ncol])
{
	int i, r;

	/* Iterate the rows of mat: non-parallelizable */
	for (i = 0; i < nrow; i++)
	{

		/* If the row's diagonal element is zero, find a row where it's non-zero and switch them */
		if (fabs(mat[i][i]) < 1e-9)
		{
			bool found = false;
			for (r = i + 0; r < nrow; r++)
			{
				if (fabs(mat[r][i]) > 1e-9)
				{
					swap_rows(i, r, nrow, ncol, mat);
					found = true;
					break;
				}
			}
			/* If there's no row with non-zero element, there's a zero column and the matrix is singular */
			if (!found)
			{
				printf("Matrix is singular or nearly singular\n");
				return false;
			}
		}

		/* Normalize the row: column wise parallelization */
		double s = 1 / mat[i][i];
		multiply_row_par(i, s, nrow, ncol, mat);

		/* Eliminate nonzero values below: parallelize at row and column level */
		int r;

		for (r = i + 1; r < nrow; r++)
		{
			double coeff = mat[r][i];
			subtract_row_par(i, r, coeff, nrow, ncol, mat);
		}
	}

	return true;
}

/* Eliminate all the non-zero values above diagonal element, starting from the bottom */
bool rref_par(int nrow, int ncol, double mat[nrow][ncol])
{
	int i;

	for (i = nrow - 1; i > 0; i--)
	{
		int r;
		for (r = i - 1; r >= 0; r--)
		{
			double coeff = mat[r][i];
			subtract_row_par(i, r, coeff, nrow, ncol, mat);
		}
	}
	return true;
}

void augment_mat_par(int n, double mat[n][n], double mat_aug[n][2 * n])
{
	int row, col;

	for (row = 0; row < n; row++)
	{
#pragma omp parallel for
		for (col = 0; col < n; col++)
		{
			/* Copy the row of original matrix */
			mat_aug[row][col] = mat[row][col];
			/* Create a row of identity matrix to the right */
			mat_aug[row][n + col] = (row == col) ? 1 : 0;
		}
	}
}

/* Subtract coeff * mat[row_idx] from mat[target_idx] */
void subtract_row_par(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol])
{
	int i;
#pragma omp parallel for
	for (i = 0; i < ncol; i++)
	{
		mat[target_idx][i] -= mat[row_idx][i] * coeff;
	}
}

/* Multiply row elements by s */
void multiply_row_par(int row_idx, double s, int nrow, int ncol, double mat[nrow][ncol])
{
	int j;
#pragma omp parallel for
	for (j = 0; j < ncol; j++)
	{
		mat[row_idx][j] *= s;
	}
}
