#include "matrix_inversion.h"
#include "matrix_inversion_parallel.h"
#include "helpers/common.h"
#include "helpers/file_reader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>

/* Helper function declarations */
double **allocate_and_read_matrix(const char *filepath, int *nrow, int *ncol);
void process_parallel_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol]);
bool invert_matrix_from_file(const char *filepath);

/* Main function to perform matrix inversion */
int main(int argc, char *argv[])
{
    printf("\n********** Matrix Inverse **********\n\n");

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

    if (!invert_matrix_from_file(filepath))
    {
        fprintf(stderr, "Failed to process file: %s\n", filepath);
        return 1;
    }

    return 0;
}

/* Function to read and invert a matrix from a file */
bool invert_matrix_from_file(const char *filepath)
{
    int nrow, ncol;
    double **mat = allocate_and_read_matrix(filepath, &nrow, &ncol);

    if (!mat)
    {
        return false;
    }

    double mat_inv_parallel[nrow][ncol];
    bool result = false;

    process_serial_inversion(nrow, ncol, mat, mat_inv_parallel);

    // printf("\n********** Inverted Matrix Start **********\n");
    // print_mat(nrow, ncol, mat_inv_parallel);
    // printf("\n********** Inverted Matrix End **********\n");

    free_matrix(mat, nrow);
    return true;
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

void process_serial_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol])
{
    double mat_cp[nrow][ncol];
    copy_matrix(nrow, ncol, mat, mat_cp);

    benchmark_matrix_inversion(nrow, ncol, mat_cp, mat_inv);
}
