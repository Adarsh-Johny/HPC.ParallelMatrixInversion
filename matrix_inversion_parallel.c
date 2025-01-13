/**
 * @file matrix_inversion_parallel.c
 * @brief Implements parallel matrix inversion using OpenMPI library
 *
 * Invert a given square matrix using gaussian elimination and benchmark the result.
 * The helper functions for augmenting the matrix, multiplying and subtracting 
 * rows and extracting the inverse are parallelized.
 * */

#include "matrix_inversion_parallel.h"
#include "helpers/common.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

/* Function for benchmarking the inversion */
void benchmark_matrix_inversion_parallel(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol])
{
	struct timeval start, end;

	// Synchronize threads before timing
#pragma omp parallel
	{
#pragma omp barrier
	}

	// Start timing
	gettimeofday(&start, NULL);

	// Perform matrix inversion
	if (!invert_matrix_par(nrow, ncol, mat, mat_inv))
	{
		printf("Matrix inversion failed during benchmarking.\n");
		return;
	}

	// Synchronize threads after execution
#pragma omp parallel
	{
#pragma omp barrier
	}

	// End timing
	gettimeofday(&end, NULL);

	// Calculate elapsed time in milliseconds
	double elapsed_time = (end.tv_sec - start.tv_sec) * 1000.0;
	elapsed_time += (end.tv_usec - start.tv_usec) / 1000.0;

	printf("Matrix inversion (Parallel) completed in %.3f ms for %dx%d matrix.\n", elapsed_time, nrow, ncol);
	// printf("Inverted Matrix:\n");
	// print_mat(nrow, ncol, mat_inv);
}

/* Invert the matrix and return the inverse */
bool invert_matrix_par(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol])
{
	int n = nrow;

	double mat_aug[n][2 * n];
	augment_mat_par(n, mat, mat_aug);

	if (!gaussian_elimination_par(n, 2 * n, mat_aug))
	{
		printf("Gaussian Elimination failed.\n");
		return false;
	}

	if (!rref_par(n, 2 * n, mat_aug))
	{
		printf("Reduced Row Echelon Form transformation failed.\n");
		return false;
	}

	extract_inverse_par(n, 2 * n, mat_aug, mat_inv);
	return true;
}

/* Extract the inverse. Input matrix is a n x 2n matrix where the right n x n matrix is the inverse */
void extract_inverse_par(int nrow, int ncol, double mat_aug[nrow][ncol], double mat_inv[nrow][nrow])
{
#pragma omp parallel for collapse(2)
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < nrow; j++)
		{
			mat_inv[i][j] = mat_aug[i][nrow + j];
		}
	}
}

/* Implementation of the gaussian elimination step. The left n x n matrix of the input matrix
 * will be a upper triangular matrix with 1s on the diagonal after this step
 */
bool gaussian_elimination_par(int nrow, int ncol, double mat[nrow][ncol])
{
	for (int i = 0; i < nrow; i++)
	{
		// Ensure the pivot element is non-zero
		if (fabs(mat[i][i]) < 1e-9)
		{
			bool swapped = false;
			for (int r = i + 1; r < nrow; r++)
			{
				if (fabs(mat[r][i]) > 1e-9)
				{
					swap_rows(i, r, nrow, ncol, mat);
					swapped = true;
					break;
				}
			}
			if (!swapped)
			{
				printf("Matrix is singular or nearly singular.\n");
				return false;
			}
		}

		// Normalize the pivot row
		double scale = 1.0 / mat[i][i];
		multiply_row_par(i, scale, nrow, ncol, mat);

// Eliminate rows below the pivot
#pragma omp parallel for schedule(dynamic)
		for (int r = i + 1; r < nrow; r++)
		{
			// printf("Thread: gaussian_elimination_par %d/%d - %d\n", omp_get_max_threads(), omp_get_thread_num(), omp_get_num_procs());
			double coeff = mat[r][i];
			subtract_row_par(i, r, coeff, nrow, ncol, mat);
		}
	}
	return true;
}

/* Turn the left n x n matrix of the input matrix into reduced row echelon form. */
bool rref_par(int nrow, int ncol, double mat[nrow][ncol])
{
	for (int i = nrow - 1; i > 0; i--)
	{
#pragma omp parallel for
		for (int r = i - 1; r >= 0; r--)
		{
			// printf("Thread rref_par: %d/%d - %d\n", omp_get_thread_num(),omp_get_max_threads(), omp_get_num_procs());

			double coeff = mat[r][i];
			subtract_row_par(i, r, coeff, nrow, ncol, mat);
		}
	}
	return true;
}

/* Create the augmented n x 2n matrix where the original input matrix is on the left with a
 * n x n identity matrix added to the right */
void augment_mat_par(int n, double mat[n][n], double mat_aug[n][2 * n])
{
#pragma omp parallel for collapse(2)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2 * n; j++)
		{
			mat_aug[i][j] = (j < n) ? mat[i][j] : (i == (j - n));
		}
	}
}

/* Subtract the row corresponding to row_idx from row with target_idx, multiplied with coefficient
 * coeff. Only the row in target_idx is altered */
void subtract_row_par(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol])
{
#pragma omp parallel for
	for (int i = 0; i < ncol; i++)
	{
		mat[target_idx][i] -= mat[row_idx][i] * coeff;
	}
}

/* Multiply the row at the given index row_idx with scale scale */
void multiply_row_par(int row_idx, double scale, int nrow, int ncol, double mat[nrow][ncol])
{
#pragma omp parallel for
	for (int j = 0; j < ncol; j++)
	{
		mat[row_idx][j] *= scale;
	}
}
