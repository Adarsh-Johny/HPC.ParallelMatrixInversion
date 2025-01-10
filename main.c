#include "matrix_inversion.h"
#include "matrix_inversion_parallel.h"
#include "helpers/common.h"
#include "helpers/file_reader.h"

#include <stdio.h>
#include <stdlib.h>  /* malloc */
#include <math.h>    /* fabs */
#include <stdbool.h> /* bool */

/* Helper function declarations */
double **allocate_and_read_matrix(const char *dir, const char *fname, int *nrow, int *ncol);
void process_serial_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result);
void process_parallel_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result);

bool invert_matrix_from_file(const char *dir, const char *fname)
{
    int nrow, ncol;
    double **mat = allocate_and_read_matrix(dir, fname, &nrow, &ncol);
    if (!mat)
    {
        return false;
    }

    double mat_inv_serial[nrow][ncol], mat_inv_parallel[nrow][ncol];
    bool serial_result = false, parallel_result = false;

    process_serial_inversion(nrow, ncol, mat, mat_inv_serial, &serial_result);
    process_parallel_inversion(nrow, ncol, mat, mat_inv_parallel, &parallel_result);

    if (serial_result && parallel_result)
    {
        compare_inversions(fname, nrow, ncol, mat_inv_serial, mat_inv_parallel);
    }
    else
    {
        if (!serial_result)
            printf("Serial inversion failed for %s\n", fname);
        if (!parallel_result)
            printf("Parallel inversion failed for %s\n", fname);
    }

    free_matrix(mat, nrow);
    printf("\n\n");
    return true;
}

double **allocate_and_read_matrix(const char *dir, const char *fname, int *nrow, int *ncol)
{
    double **mat;
    if (!read_matrix_from_file(dir, fname, nrow, ncol, &mat))
    {
        fprintf(stderr, "Failed to read matrix from file %s\n", fname);
        return NULL;
    }
    return mat;
}

void process_serial_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result)
{
    double mat_cp[nrow][ncol];
    copy_matrix(nrow, ncol, mat, mat_cp);
    *result = invert_matrix(nrow, ncol, mat_cp, mat_inv);
}

void process_parallel_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result)
{
    double mat_cp[nrow][ncol];
    copy_matrix(nrow, ncol, mat, mat_cp);
    *result = invert_matrix_par(nrow, ncol, mat_cp, mat_inv);
}
