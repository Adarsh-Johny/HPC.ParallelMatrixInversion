#include "matrix_inversion.h"
#include "helpers/common.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h> /* gettimeofday */
#include <math.h>     /* fabs */

#define EPSILON 1e-9
#define DEBUG 1
#if DEBUG
#define DBG_PRINT(...) printf(__VA_ARGS__)
#else
#define DBG_PRINT(...)
#endif

#include <sys/time.h> /* gettimeofday */

void benchmark_matrix_inversion(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol])
{
    struct timeval start, end;

    /* Start timing */
    gettimeofday(&start, NULL);

    /* Call the matrix inversion function */
    if (!invert_matrix(nrow, ncol, mat, mat_inv))
    {
        printf("Matrix inversion failed during benchmarking.\n");
        return;
    }

    /* End timing */
    gettimeofday(&end, NULL);

    /* Calculate elapsed time in milliseconds */
    double elapsed_time = (end.tv_sec - start.tv_sec) * 1000.0; /* Seconds to milliseconds */
    elapsed_time += (end.tv_usec - start.tv_usec) / 1000.0;     /* Microseconds to milliseconds */

    printf("Matrix inversion (Serial) completed in %.3f ms for %dx%d matrix.\n", elapsed_time, nrow, ncol);
}

bool invert_matrix(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol])
{
    if (nrow != ncol)
    {
        printf("Matrix inversion only supported for square matrices.\n");
        return false;
    }

    int n = nrow;

    /* Augment identity */
    double mat_aug[n][2 * n];
    augment_mat_ser(n, mat, mat_aug);

    /* Forward elimination */
    if (!gaussian_elimination(n, 2 * n, mat_aug))
    {
        printf("Gaussian elimination failed.\n");
        return false;
    }
    // DBG_PRINT("Matrix after Gaussian elimination:\n");
    // print_mat(n, 2 * n, mat_aug);

    /* Backward elimination */
    if (!rref(n, 2 * n, mat_aug))
    {
        printf("RREF failed.\n");
        return false;
    }
    // DBG_PRINT("Matrix after RREF:\n");
    // print_mat(n, 2 * n, mat_aug);

    /* Extract inverse */
    extract_inverse_ser(n, 2 * n, mat_aug, mat_inv);

    // DBG_PRINT("Extracted Inverse Matrix:\n");
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         printf("%8.9f ", mat_inv[i][j]);
    //     }
    //     printf("\n");
    // }

    return true;
}

void multiply_row_ser(int row_idx, double s, int nrow, int ncol, double mat[nrow][ncol])
{
    if (row_idx < 0 || row_idx >= nrow)
    {
        printf("Error: Invalid row index %d for matrix %d x %d\n", row_idx, nrow, ncol);
        return;
    }

    for (int j = 0; j < ncol; j++)
    {
        mat[row_idx][j] *= s;
    }
}

void augment_mat_ser(int n, const double mat[n][n], double mat_aug[n][2 * n])
{
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            mat_aug[row][col] = mat[row][col];
            mat_aug[row][n + col] = (row == col) ? 1 : 0;
        }
    }
}

bool gaussian_elimination(int nrow, int ncol, double mat[nrow][ncol])
{
    for (int i = 0; i < nrow; i++)
    {
        if (fabs(mat[i][i]) < EPSILON)
        {
            bool found = false;
            for (int r = i + 1; r < nrow; r++)
            {
                if (fabs(mat[r][i]) > EPSILON)
                {
                    swap_rows(i, r, nrow, ncol, mat);
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                printf("Matrix is singular at row %d, column %d\n", i, i);
                return false;
            }
        }

        double s = 1.0 / mat[i][i];
        multiply_row_ser(i, s, nrow, ncol, mat);

        for (int r = i + 1; r < nrow; r++)
        {
            double coeff = mat[r][i];
            subtract_row_ser(i, r, coeff, nrow, ncol, mat);
        }
    }
    return true;
}

void subtract_row_ser(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol])
{
    if (row_idx < 0 || row_idx >= nrow)
    {
        printf("Error: Invalid source row index %d for matrix %d x %d\n", row_idx, nrow, ncol);
        return;
    }

    if (target_idx < 0 || target_idx >= nrow)
    {
        printf("Error: Invalid target row index %d for matrix %d x %d\n", target_idx, nrow, ncol);
        return;
    }

    for (int j = 0; j < ncol; j++)
    {
        mat[target_idx][j] -= coeff * mat[row_idx][j];
    }
}

bool rref(int nrow, int ncol, double mat[nrow][ncol])
{
    for (int i = nrow - 1; i > 0; i--)
    {
        for (int r = i - 1; r >= 0; r--)
        {
            double coeff = mat[r][i];
            subtract_row_ser(i, r, coeff, nrow, ncol, mat);
        }
    }
    return true;
}

void extract_inverse_ser(int nrow, int ncol, const double mat_aug[nrow][ncol], double mat_inv[nrow][nrow])
{
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < nrow; j++)
        {
            mat_inv[i][j] = mat_aug[i][nrow + j];
        }
    }
}
