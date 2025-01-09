#include "matrix_inversion.h"
#include "matrix_inversion_parallel.h"
#include "common.h"

#include <stdio.h>
#include <stdlib.h> /* malloc */
#include <math.h>   /* fabs */
#include <dirent.h> /* DIR, readdir */
#include <string.h> /* strlen */
#include <sys/types.h>
#include <sys/stat.h>

/* Run test matrices */
bool invert_matrix_from_file(const char *dir, const char *fname);
bool check_inverse(int nrow, int ncol, double m1[nrow][ncol], double m2[nrow][ncol]);

int main(int argc, char *argv[])
{
    /* int n = 3; */ /* Size of the matrix */

    /* Get and print input matrix */
    /*
       double mat[3][3] = { {0, 3, 2}, {1, 0, 2}, {0, 0, 1} };
       int nrow = sizeof(mat) / sizeof(mat[0]);
       int ncol = sizeof(mat[0]) / sizeof(mat[0][0]);

       print_mat(nrow, ncol, mat);
       */

    // FILE *fp = fopen("HPC.ParallelMatrixInversion/test_matrices/mat_3x3_i_1.txt", "r");
    // /* Check that the file exists */
    // if (!fp)
    // {
    //     perror("fopen");
    //     return 1;
    // }

    const char *directory_path = "HPC.ParallelMatrixInversion/test_matrices/";
    struct dirent *entry;
    DIR *dir = opendir(directory_path);

    if (!dir)
    {
        perror("opendir");
        return 1;
    }

    printf("Opened dir %s\n", directory_path);

    while ((entry = readdir(dir)) != NULL)
    {
        // Construct the full path to the file
        char full_path[1024];
        snprintf(full_path, sizeof(full_path), "%s/%s", directory_path, entry->d_name);

        // Use stat() to check if it's a regular file
        struct stat file_stat;
        if (stat(full_path, &file_stat) == 0 && S_ISREG(file_stat.st_mode))
        {
            invert_matrix_from_file(directory_path, entry->d_name);
        }
    }

    closedir(dir);

    return 0;
}

/* A method to run test matrices on invert_matrix */
bool read_matrix_from_file(const char *dir, const char *fname, int *nrow, int *ncol, double ***mat)
{
    /* Create the target file path */
    char *full_path = malloc(strlen(dir) + strlen(fname) + 2); // +2 for '/' and '\0'

    if (!full_path)
    {
        perror("malloc (concatenate file path)");
        return false;
    }

    strcpy(full_path, dir);
    strcat(full_path, "/");
    strcat(full_path, fname);

    /* Extract dimensions from the filename */
    printf("Reading matrix from file %s\n", fname);
    int index;

    if (sscanf(fname, "matrix_%dx%d_%02d.txt", nrow, ncol, &index) != 3)
    {
        printf("Invalid filename format: %s\n", fname);
        free(full_path);
        return false;
    }

    printf("Reading %dx%d matrix from %s\n", *nrow, *ncol, full_path);

    /* Open the file */
    FILE *fp = fopen(full_path, "r");
    if (!fp)
    {
        perror("Error opening file");
        free(full_path);
        return false;
    }

    /* Allocate memory for the matrix */
    *mat = (double **)malloc(*nrow * sizeof(double *));
    for (int i = 0; i < *nrow; ++i)
    {
        (*mat)[i] = (double *)malloc(*ncol * sizeof(double));
    }

    /* Scan the file and read the matrix */
    for (int i = 0; i < *nrow; ++i)
    {
        for (int j = 0; j < *ncol; ++j)
        {
            if (fscanf(fp, "%lf", &(*mat)[i][j]) != 1)
            {
                printf("Error reading matrix value at [%d][%d] in file: %s\n", i, j, fname);
                fclose(fp);
                free(full_path);
                return false;
            }
        }
    }

    fclose(fp);
    free(full_path);
    return true;
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

bool invert_matrix_and_validate(int nrow, int ncol, double source[nrow][ncol], double inverse[nrow][ncol])
{
    bool res = invert_matrix(nrow, ncol, source, inverse);

    if (res)
    {
        if (check_inverse(nrow, ncol, source, inverse))
        {
            printf("Found correct inverse\n");
        }
        else
        {
            printf("ERROR: Didn't find inverse for matrix\n");
        }
    }
    else
    {
        printf("ERROR: Failed to invert matrix\n");
    }

    return res;
}

bool invert_matrix_from_file(const char *dir, const char *fname)
{
    int nrow, ncol;
    double **mat;

    if (!read_matrix_from_file(dir, fname, &nrow, &ncol, &mat))
    {
        return false;
    }

    /* Create a copy of the input matrix before inverting it */
    double mat_cp[nrow][ncol];
    copy_matrix(nrow, ncol, mat, mat_cp);

    /* Invert the matrix */
    double mat_inv[nrow][ncol];
    invert_matrix_and_validate(nrow, ncol, mat_cp, mat_inv);

    /* Free allocated memory */
    for (int i = 0; i < nrow; ++i)
    {
        free(mat[i]);
    }
    free(mat);

    printf("\n\n");
    return true;
}

bool check_inverse(int nrow, int ncol, double m1[nrow][ncol], double m2[nrow][ncol])
{
    bool success = true;
    double mat_res[nrow][ncol];
    int i, j, k;

    for (i = 0; i < nrow; i++)
    {
        for (j = 0; j < ncol; j++)
        {
            mat_res[i][j] = 0;

            for (k = 0; k < nrow; k++)
            {
                mat_res[i][j] += m1[i][k] * m2[k][j];
            }

            /* Check if the result is identity (will fail if there are nans or +-inf) */
            if (i == j && fabs(1 - mat_res[i][j]) > 1e-9)
            {
                success = false;
            }

            if (i != j && fabs(mat_res[i][j]) > 1e-9)
            {
                success = false;
            }
        }
    }

    return success;
}
