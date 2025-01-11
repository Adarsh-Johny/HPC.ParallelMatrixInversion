#include "matrix_inversion.h"
#include "matrix_inversion_parallel.h"
#include "helpers/common.h"
#include "helpers/file_reader.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>

/* Helper function declarations */
double **allocate_and_read_matrix(const char *dir, const char *fname, int *nrow, int *ncol);
void process_parallel_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result);
bool invert_matrix_from_file(const char *dir, const char *fname);
void process_serial_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result);

/* Main function to perform matrix inversion */
int main(int argc, char *argv[])
{
    const char *directory_path = "HPC.ParallelMatrixInversion/performance_test_matrices/";
    DIR *dir = opendir(directory_path);

    if (!dir)
    {
        perror("opendir");
        return 1;
    }

    printf("Opened directory: %s\n", directory_path);

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL)
    {
        if (entry->d_name[0] == '.') // Skip hidden files and "." or ".."
            continue;

        // Check if the file matches the expected format
        if (strstr(entry->d_name, "matrix_") == NULL)
            continue;

        // Construct the full path and process
        char full_path[1024];
        snprintf(full_path, sizeof(full_path), "%s/%s", directory_path, entry->d_name);

        struct stat file_stat;
        if (stat(full_path, &file_stat) == 0 && S_ISREG(file_stat.st_mode))
        {
            invert_matrix_from_file(directory_path, entry->d_name);
        }
    }

    closedir(dir);
    return 0;
}

/* Function to read and invert a matrix from a file */
bool invert_matrix_from_file(const char *dir, const char *fname)
{
    fprintf("Called file: %s", *fname);

    int nrow, ncol;
    double **mat = allocate_and_read_matrix(dir, fname, &nrow, &ncol);

    if (!mat)
    {
        return false;
    }

    // printf("\n********** Original Matrix (%s) **********\n", fname);
    // print_mat(nrow, ncol, mat);

    double mat_inv_parallel[nrow][ncol];
    bool result = false;

    process_parallel_inversion(nrow, ncol, mat, mat_inv_parallel, &result);
    // process_serial_inversion(nrow, ncol, mat, mat_inv_parallel, &result);

    // if (result)
    // {
    // printf("\n********** Inverted Matrix (Parallel) **********\n");
    // print_mat(nrow, ncol, mat_inv_parallel);
    // }
    // else
    // {
    //     printf("Matrix inversion (Parallel) failed for file: %s\n", fname);
    // }

    free_matrix(mat, nrow);
    return true;
}

/* Helper function to allocate and read a matrix */
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

/* Process parallel matrix inversion */
void process_parallel_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result)
{
    double mat_cp[nrow][ncol];
    copy_matrix(nrow, ncol, mat, mat_cp);

    // Benchmark and invert matrix
    benchmark_matrix_inversion_parallel(nrow, ncol, mat_cp, mat_inv);
    // *result = invert_matrix_par(nrow, ncol, mat_cp, mat_inv);
}

/* Process parallel matrix inversion */
void process_serial_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result)
{
    double mat_cp[nrow][ncol];
    copy_matrix(nrow, ncol, mat, mat_cp);

    // Benchmark and invert matrix
    benchmark_matrix_inversion(nrow, ncol, mat_cp, mat_inv);
    // *result = invert_matrix_par(nrow, ncol, mat_cp, mat_inv);
}
