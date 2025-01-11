#include <stdio.h>     /* printf (if needed elsewhere) */
#include <stdlib.h>    /* malloc (if used elsewhere) */
#include <unistd.h>    /* getcwd */
#include <math.h>      /* fabs */
#include <dirent.h>    /* DIR, readdir (if directory traversal is needed elsewhere) */
#include <string.h>    /* strlen (if string operations are needed) */
#include <sys/types.h> /* stat */
#include <sys/stat.h>  /* stat */
#include <stdbool.h>   /* bool, true, false */

/* Swap rows r1 and r2 */
void swap_rows(int r1, int r2, int nrow, int ncol, double mat[nrow][ncol])
{
    double temp;
    int i;
    for (i = 0; i < ncol; i++)
    {
        temp = mat[r1][i];
        mat[r1][i] = mat[r2][i];
        mat[r2][i] = temp;
    }
}

void print_working_dir()
{
    /* Check the working directory */
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
    {
        printf("Current working directory: %s\n", cwd);
    }
    else
    {
        perror("getcwd");
    }
}

void print_mat(int nrow, int ncol, double mat[nrow][ncol])
{
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            printf("%8.12f ", mat[i][j]);
        }
        printf("\n");
    }
}

void copy_mat(int nrow, int ncol, double src[nrow][ncol], double dest[nrow][ncol])
{
    int i, j;
    for (i = 0; i < nrow; i++)
    {
        for (j = 0; j < ncol; j++)
        {
            dest[i][j] = src[i][j];
        }
    }
}

void copy_matrix(int nrow, int ncol, double **source, double dest[nrow][ncol])
{
    for (int i = 0; i < nrow; ++i)
    {
        for (int j = 0; j < ncol; ++j)
        {
            dest[i][j] = source[i][j];
        }
    }
}

bool compare_matrices(int nrow, int ncol, double mat1[nrow][ncol], double mat2[nrow][ncol])
{
    const double tolerance = 1e-6; // Adjust as needed
    for (int i = 0; i < nrow; ++i)
    {
        for (int j = 0; j < ncol; ++j)
        {
            if (fabs(mat1[i][j] - mat2[i][j]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}
bool check_inverse(int nrow, int ncol, double m1[nrow][ncol], double m2[nrow][ncol])
{
    if (nrow != ncol)
    {
        printf("Error: Matrices must be square to check for inverses.\n");
        return false;
    }

    bool success = true;
    double mat_res[nrow][ncol];
    double tolerance = 1e-9;
    int i, j, k;

    // Initialize mat_res to zero and compute the product of m1 and m2
    for (i = 0; i < nrow; i++)
    {
        for (j = 0; j < ncol; j++)
        {
            mat_res[i][j] = 0;

            for (k = 0; k < nrow; k++)
            {
                mat_res[i][j] += m1[i][k] * m2[k][j];
            }

            // Check identity matrix conditions
            if (i == j && fabs(1.0 - mat_res[i][j]) > tolerance)
            {
                printf("Diagonal element at (%d, %d) is not 1. Value: %f\n", i, j, mat_res[i][j]);
                success = false;
            }
            else if (i != j && fabs(mat_res[i][j]) > tolerance)
            {
                printf("Off-diagonal element at (%d, %d) is not 0. Value: %f\n", i, j, mat_res[i][j]);
                success = false;
            }
        }
    }

    printf("####check_inverse: original matrix and its inverse multiplication:\n");
    print_mat(nrow, ncol, mat_res);

    return success;
}

void compare_inversions(const char *fname, int nrow, int ncol, double mat_inv_serial[nrow][ncol], double mat_inv_parallel[nrow][ncol])
{
    if (compare_matrices(nrow, ncol, mat_inv_serial, mat_inv_parallel))
    {
        printf("Both implementations produced the same results for %s\n", fname);
    }
    else
    {
        printf("Mismatch between serial and parallel results for %s\n", fname);
    }
}

void free_matrix(double **mat, int nrow)
{
    for (int i = 0; i < nrow; ++i)
    {
        free(mat[i]);
    }
    free(mat);
}
