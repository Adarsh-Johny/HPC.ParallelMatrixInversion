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
    double **augmented = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        augmented[i] = malloc(2 * n * sizeof(double));
    }

    // Initialize augmented matrix [mat | I]
    if (rank == 0)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                augmented[i][j] = mat[i][j];
                augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    // Broadcast augmented matrix to all processes
    for (int i = 0; i < n; i++)
    {
        MPI_Bcast(augmented[i], 2 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Gaussian elimination
    for (int k = 0; k < n; k++)
    {
        if (fabs(augmented[k][k]) < EPSILON)
        {
            if (rank == 0)
            {
                fprintf(stderr, "Matrix is singular or nearly singular.\n");
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Scale pivot row
        if (rank == k % size)
        {
            double pivot = augmented[k][k];
            for (int j = 0; j < 2 * n; j++)
            {
                augmented[k][j] /= pivot;
            }
        }

        MPI_Bcast(augmented[k], 2 * n, MPI_DOUBLE, k % size, MPI_COMM_WORLD);

        // Eliminate other rows
        for (int i = 0; i < n; i++)
        {
            if (i != k)
            {
                double factor = augmented[i][k];
                for (int j = 0; j < 2 * n; j++)
                {
                    augmented[i][j] -= factor * augmented[k][j];
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
                mat_inv_parallel[i][j] = augmented[i][j + n];
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        free(augmented[i]);
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        double elapsed_time = (end_time - start_time) * 1000.0; // Convert seconds to milliseconds

        printf("Matrix inversion (Parallel) completed in %.3f ms for %dx%d matrix.\n", elapsed_time, nrow, ncol);

        // printf("Inverted Matrix:\n");
        // print_mat_pointer(nrow, ncol, mat_inv_parallel);
    }
}