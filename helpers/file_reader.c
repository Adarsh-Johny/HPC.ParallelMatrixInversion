#include "file_reader.h"
#include <stdio.h>   /* printf, perror, FILE, fopen, fscanf */
#include <stdlib.h>  /* malloc, free */
#include <string.h>  /* strlen, strcpy, strcat */
#include <stdbool.h> /* bool, true, false */

/* Helper function to free the matrix */
void free_matrix_file_reader(double **mat, int nrow)
{
    if (mat)
    {
        for (int i = 0; i < nrow; ++i)
        {
            free(mat[i]);
        }
        free(mat);
    }
}

/* A method to read a matrix from a file */
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
    if (!*mat)
    {
        perror("malloc (matrix rows)");
        fclose(fp);
        free(full_path);
        return false;
    }

    for (int i = 0; i < *nrow; ++i)
    {
        (*mat)[i] = (double *)malloc(*ncol * sizeof(double));
        if (!(*mat)[i])
        {
            perror("malloc (matrix columns)");
            free_matrix_file_reader(*mat, i); // Free already allocated rows
            fclose(fp);
            free(full_path);
            return false;
        }
    }

    /* Scan the file and read the matrix */
    for (int i = 0; i < *nrow; ++i)
    {
        for (int j = 0; j < *ncol; ++j)
        {
            if (fscanf(fp, "%lf", &(*mat)[i][j]) != 1)
            {
                printf("Error reading matrix value at [%d][%d] in file: %s\n", i, j, fname);
                free_matrix_file_reader(*mat, *nrow);
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
