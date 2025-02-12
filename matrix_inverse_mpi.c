/*
 * @file matrix_inverse_mpi.c
 * @brief Implements a parallel matrix inversion using MPI
 *
 * Inverse the given matrix with inverse_matrix_mpi, called from the benchmarking function
 * benchmark_inversion to measure the performance of the parallel implementation.
 * Functions expect that the input is a square matrix.
 */

#include "matrix_inverse_mpi.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define EPSILON 1e-10

void inverse_matrix_mpi(double **mat, int nrow, int ncol, double **mat_inv_parallel)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (nrow != ncol)
    {
        if (rank == 0)
        {
            fprintf(stderr, "Matrix must be square for inversion.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int n = nrow;
    double *augmented = malloc(n * 2 * n * sizeof(double)); // Contiguous memory

    // Initialize augmented matrix [mat | I]
    if (rank == 0)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                augmented[i * 2 * n + j] = mat[i][j];
                augmented[i * 2 * n + j + n] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    // Broadcast augmented matrix to all processes
    MPI_Bcast(augmented, n * 2 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Gaussian elimination
    for (int k = 0; k < n; k++)
    {
        int pivot_rank = k % size;
        if (rank == pivot_rank)
        {
            if (fabs(augmented[k * 2 * n + k]) < EPSILON)
            {
                fprintf(stderr, "Matrix is singular or nearly singular.\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // Scale pivot row
            double pivot = augmented[k * 2 * n + k];
            for (int j = 0; j < 2 * n; j++)
            {
                augmented[k * 2 * n + j] /= pivot;
            }
        }

        MPI_Bcast(&augmented[k * 2 * n], 2 * n, MPI_DOUBLE, pivot_rank, MPI_COMM_WORLD);

        // Eliminate other rows
        for (int i = 0; i < n; i++)
        {
            if (i != k)
            {
                double factor = augmented[i * 2 * n + k];
                for (int j = 0; j < 2 * n; j++)
                {
                    augmented[i * 2 * n + j] -= factor * augmented[k * 2 * n + j];
                }
            }
        }
    }

    // Extract inverse matrix
    if (rank == 0)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                mat_inv_parallel[i][j] = augmented[i * 2 * n + j + n];
            }
        }
    }

    free(augmented);
}

void benchmark_inversion(double **mat, int nrow, int ncol)
{
    double **mat_inv_parallel = malloc(nrow * sizeof(double *));
    for (int i = 0; i < nrow; i++)
    {
        mat_inv_parallel[i] = malloc(ncol * sizeof(double));
    }

    double start_time = MPI_Wtime();

    inverse_matrix_mpi(mat, nrow, ncol, mat_inv_parallel);

    double end_time = MPI_Wtime();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Determines the Rank of process

    if (rank == 0)
    {
        double elapsed_time = (end_time - start_time) * 1000.0; // Convert seconds to milliseconds

        printf("Matrix inversion (Parallel) completed in %.3f ms for %dx%d matrix.\n", elapsed_time, nrow, ncol);

        // printf("Inverted Matrix:\n");
        // print_mat_pointer(nrow, ncol, mat_inv_parallel);
    }
}
