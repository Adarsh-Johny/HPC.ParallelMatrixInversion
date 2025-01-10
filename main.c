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
#include <string.h> /* strlen, strcpy, strcat */

/* Helper function declarations */
double **allocate_and_read_matrix(const char *dir, const char *fname, int *nrow, int *ncol);
void process_serial_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result);
void process_parallel_inversion(int nrow, int ncol, double **mat, double mat_inv[nrow][ncol], bool *result);
bool invert_matrix_from_file(const char *dir, const char *fname);

bool invert_test_mat_from_file(const char* dir, const char* fname, bool test_parallel);
void run_tests();

int main(int argc, char *argv[]) {

    run_tests();
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

void run_tests() {

    const char *directory_path = "HPC.ParallelMatrixInversion/test_matrices/";
    DIR *dir = opendir(directory_path);
    if (!dir) {
	perror("opendir");
	printf("Couldn't open directory %s\n", directory_path);
	return;
    }   

    printf("Opened dir %s\n", directory_path);


    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) 
    {
	/* Skip directories (there shouldn't be any) */
	/* if (entry->d_type == DT_REG) */
	char full_path[1024];
        snprintf(full_path, sizeof(full_path), "%s/%s", directory_path, entry->d_name);

	struct stat file_stat;
        if (stat(full_path, &file_stat) == 0 && S_ISREG(file_stat.st_mode))
	{ 
	    if (!invert_test_mat_from_file(directory_path, entry->d_name, true))
	    {
		printf("Running tests for file %s failed (parallel) \n", entry->d_name);
	    }
	    /*
	    if (!invert_test_mat_from_file(directory_path, entry->d_name, false))
	    {		
		printf("Running tests for file %s failed (sequential) \n", entry->d_name);
	    }
	    */
	}   
    }   
    closedir(dir);;
}   


/* A method to run test matrices on invert_matrix */
bool invert_test_mat_from_file(const char* dir, const char* fname, bool test_parallel) {
    /* Create the target file path  */
    char* full_path = malloc(strlen(dir) + strlen(fname) + 1);

    if (!full_path) {
	perror("malloc (concatenate file path)");
	return false;
    }   
    strcpy(full_path, dir);
    strcat(full_path, fname);

    /* Extract dimensions and type of matrix (invertible or singular) from the filename */
    printf("Reading matrix from file %s\n", fname);
    int nrow, ncol;
    char kind;

    if (sscanf(fname, "matrix_%dx%d_%c", &nrow, &ncol, &kind) != 3) {
	printf("Invalid filename format: %s\n", fname);
	return false;
    }

    printf("Reading %dx%d matrix from %s ", nrow, ncol, full_path);
    printf("(%s)\n", (kind == 'i') ? "invertible" : "singular");

    /* Open the file */
    FILE *fp = fopen(full_path, "r");
    if (!fp) {
	perror("Error opening file");
	print_working_dir();
	return false;
    }

    /* Scan the file and read the matrix */
    int i, j;
    double mat[nrow][ncol];

    for (i = 0; i < nrow; ++i) {
	for (j = 0; j < ncol; ++j) {
	    if (fscanf(fp, "%lf", &mat[i][j]) != 1) {
		printf("Error reading matrix value at [%d][%d] in file: %s\n", i, j, fname);
		fclose(fp);
		return false;
	    }
	}
    }

    fclose(fp);
    free(full_path);

    /* Create a copy of the input matrix before inverting it */
    double mat_cp[nrow][ncol];
    copy_mat(nrow, ncol, mat, mat_cp);

    /* Invert the matrix */
    double mat_inv[nrow][ncol];
    bool res;
    if (test_parallel)
    {
	res = invert_matrix_par(nrow, ncol, mat_cp, mat_inv);
    }
    else
    {
	res = invert_matrix(nrow, ncol, mat_cp, mat_inv);
    }


    if (res && kind == 'i') {

	printf("Found inverse:\n");
        print_mat(nrow, ncol, mat_inv);

	if (check_inverse(nrow, ncol, mat, mat_inv)) {
	    printf("Found correct inverse\n");
	} else {
	    printf("ERROR: Didn't find inverse for invertible matrix\n");
	}
    } else if (!res && kind == 'i') {
	printf("ERROR: Failed to invert invertible matrix\n");
    } else if (!res && kind == 's') {
	printf("Matrix is singular\n");
    }

    printf("---------------------\n\n");
    return true;
}

