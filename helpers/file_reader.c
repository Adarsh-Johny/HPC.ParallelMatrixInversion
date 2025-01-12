#include "file_reader.h"
#include <stdio.h>   /* printf, perror, FILE, fopen, fscanf */
#include <stdlib.h>  /* malloc, free */
#include <string.h>  /* strlen, strcpy, strcat */
#include <stdbool.h> /* bool, true, false */
#include <stddef.h>

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
bool read_matrix_from_file(const char *filepath, int *nrow, int *ncol, double ***mat)
{
    /* Extract dimensions from the filename */
    int index;
    const char *fname = strrchr(filepath, '/'); // Extract filename from path
    if (!fname)
        fname = filepath; // If no '/' found, the entire path is the filename
    else
        fname++; // Move past the '/'

    if (sscanf(fname, "matrix_%dx%d_%02d.txt", nrow, ncol, &index) != 3)
    {
        printf("Skipping invalid filename: %s\n", fname);
        return false;
    }

    printf("Reading %dx%d matrix from %s\n", *nrow, *ncol, filepath);

    /* Open the file */
    FILE *fp = fopen(filepath, "r");
    if (!fp)
    {
        perror("Error opening file");
        return false;
    }

    /* Allocate memory for the matrix */
    *mat = (double **)malloc(*nrow * sizeof(double *));
    if (!*mat)
    {
        perror("malloc (matrix rows)");
        fclose(fp);
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
                printf("Error reading matrix value at [%d][%d] in file: %s\n", i, j, filepath);
                free_matrix_file_reader(*mat, *nrow);
                fclose(fp);
                return false;
            }
        }
    }

    fclose(fp);
    return true;
}
