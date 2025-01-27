// main.c
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix_inverse_mpi.h"
#include "helpers/common.h"
#include "helpers/file_reader.h"
#include <string.h>
#include <stdbool.h>
#include <stddef.h>

#define EPSILON 1e-10

double **allocate_and_read_matrix(const char *filepath, int *nrow, int *ncol);
bool invert_matrix_from_file(const char *filepath);

/* Main function to perform matrix inversion */
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv); // Initializes the MPI environment and sets up communication between processes.

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s -path=<file_path>\n", argv[0]);
        return 1;
    }

    const char *filepath = NULL;

    // Parse command-line arguments
    for (int i = 1; i < argc; i++)
    {
        if (strncmp(argv[i], "-path=", 6) == 0)
        {
            filepath = argv[i] + 6; // Extract file path
        }
    }

    if (!filepath)
    {
        fprintf(stderr, "Error: File path not specified. Use -path=<file_path>.\n");
        return 1;
    }

    int nrow, ncol;
    double **mat = allocate_and_read_matrix(filepath, &nrow, &ncol);

    // printf("\n********** Matrix Original start**********\n\n");

    // print_mat_pointer(nrow, ncol, mat);

    // printf("\n********** Matrix Original end**********\n\n");

    if (!mat)
    {
        return false;
    }

    benchmark_inversion(mat, nrow, ncol);

    for (int i = 0; i < nrow; i++)
    {
        free(mat[i]);
    }
    free(mat);

    MPI_Finalize(); // Clean up all resources allocated

    return 0;
}

/* Helper function to allocate and read a matrix */
double **allocate_and_read_matrix(const char *filepath, int *nrow, int *ncol)
{
    double **mat;
    if (!read_matrix_from_file(filepath, nrow, ncol, &mat))
    {
        fprintf(stderr, "Failed to read matrix from file %s\n", filepath);
        return NULL;
    }
    return mat;
}
