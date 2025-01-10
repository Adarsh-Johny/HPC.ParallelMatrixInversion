#include "matrix_inversion.h"
#include "matrix_inversion_parallel.h"
#include "helpers/common.h"
#include "helpers/file_reader.h"

#include <stdio.h>
#include <stdlib.h>  /* malloc */
#include <math.h>    /* fabs */
#include <stdbool.h> /* bool */

#include <dirent.h> /* DIR, readdir */
#include <sys/types.h> /* struct stat */
#include <sys/stat.h> /* stat, IS_REG */
/* #include <unistd.h> */

/* Helper function declarations */
double **allocate_and_read_matrix(const char *dir, const char *fname, int *nrow, int *ncol);
void process_serial_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result);
void process_parallel_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result);
bool invert_matrix_from_file(const char *dir, const char *fname);

int main(int argc, char *argv[]) {

    const char *directory_path = "HPC.ParallelMatrixInversion/performance_test_matrices/";
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
	/* Construct the full path to the file */
	char full_path[1024];
	snprintf(full_path, sizeof(full_path), "%s/%s", directory_path, entry->d_name);

	/* Use stat() to check if it's a regular file */
	struct stat file_stat;
	if (stat(full_path, &file_stat) == 0 && S_ISREG(file_stat.st_mode))
	{
	    if (!invert_matrix_from_file(directory_path, entry->d_name))
	    {
		printf("Failed to invert matrix from %s\n", full_path);
	    }
	}
    }

    closedir(dir);

    return 0;

}

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
