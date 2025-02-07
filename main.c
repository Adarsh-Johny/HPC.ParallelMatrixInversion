#include "matrix_inversion.h"
#include "matrix_inversion_parallel.h"
#include "common.h"

#include <stdio.h>
#include <stdlib.h> /* malloc */
#include <math.h> /* fabs */
#include <dirent.h> /* DIR, readdir */
#include <string.h> /* strlen */

/* Run test matrices */
bool invert_mat_from_file(const char* dir, const char* fname);
bool check_inverse(int nrow, int ncol, double m1[nrow][ncol], double m2[nrow][ncol]);


int main(int argc, char* argv[]) {
    /* int n = 3; */ /* Size of the matrix */

    /* Get and print input matrix */
    /*
       double mat[3][3] = { {0, 3, 2}, {1, 0, 2}, {0, 0, 1} };
       int nrow = sizeof(mat) / sizeof(mat[0]);
       int ncol = sizeof(mat[0]) / sizeof(mat[0][0]);

       print_mat(nrow, ncol, mat);
       */


    // FILE *fp = fopen("HPC.ParallelMatrixInversion/test_matrices/mat_3x3_i_1.txt", "r");
    /* Check that the file exists */
    /*
    if (!fp) {
	perror("fopen");
	return 1;
    }
    */    

    const char *directory_path = "HPC.ParallelMatrixInversion/test_matrices/";
    DIR *dir = opendir(directory_path);
    if (!dir) {
	perror("opendir");
	return 1;
    }

    printf("Opened dir %s\n", directory_path);


    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
	/* Skip directories (there shouldn't be any) */
	if (entry->d_type == DT_REG) {
	    invert_mat_from_file(directory_path, entry->d_name);
	}
    }
    closedir(dir);

    return 0;
}


/* A method to test invert_matrix with specific set of test matrices 
 * The test matrices are either invertible or singular, which is determined by their
 * filename: the format is mat_nxn_c, where n is the dimension and c is either i (invertible)
 *  or s (singular)
*/
bool invert_mat_from_file(const char* dir, const char* fname) {

    /* Create the target file path  */
    char* full_path = malloc(strlen(dir) + strlen(fname) + 1);

    if (!full_path) {
	perror("malloc (concatenate file path)");
	return false;
    }

    strcpy(full_path, dir);
    strcat(full_path, fname);

    /* Extract dimensions and type of matrix (invertible or singular) from the filename */
    printf("Reading matrix from file %s ", fname);
    int nrow, ncol;
    char kind;

    if (sscanf(fname, "mat_%dx%d_%c", &nrow, &ncol, &kind) != 3) {
	printf("Invalid filename format: %s\n", fname);
	return false;
    }

    // printf("Reading %dx%d matrix from %s ", nrow, ncol, full_path);
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

    /* Close the file after successfully reading the matrix */
    fclose(fp);
    free(full_path);

    /* Create a copy of the input matrix before inverting it */
    double mat_cp[nrow][ncol];
    copy_mat(nrow, ncol, mat, mat_cp);

    /* Invert the copy of original matrix */
    double mat_inv[nrow][ncol];
    bool res = invert_matrix(nrow, ncol, mat_cp, mat_inv);

    if (res && kind == 'i') {
	if (check_inverse(nrow, ncol, mat, mat_inv)) {
	    printf("Found correct inverse for %s\n", fname);
	} else {
	    printf("ERROR: Didn't find inverse for invertible matrix at %s\n", fname);
	}
    } else if (!res && kind == 'i') {
	printf("ERROR: Failed to invert invertible matrix at %s\n", fname);
    }

    printf("\n\n");
    return true;
}

/* Function to test that the found inverse is the actual inverse. Used as a test for generated
 * matrices */
bool check_inverse(int nrow, int ncol, double m1[nrow][ncol], double m2[nrow][ncol]) {
    bool success = true;
    double mat_res[nrow][ncol];
    int i, j, k;

    for (i = 0; i < nrow; i++) {
	for (j = 0; j < ncol; j++) {
	    mat_res[i][j] = 0;

	    for (k = 0; k < nrow; k++) {
		mat_res[i][j] += m1[i][k] * m2[k][j];
	    }

	    /* Check if the result is identity (will fail if there are nans or +-inf) */
	    if (i == j && fabs(1 - mat_res[i][j]) > 1e-9) {
		success = false;
	    }

	    if (i != j && fabs(mat_res[i][j]) > 1e-9) {
		success = false;
	    }
	}
    }

    return success;
}

