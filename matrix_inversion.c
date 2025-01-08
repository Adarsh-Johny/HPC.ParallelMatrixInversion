#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* getcwd */
#include <stdbool.h>
#include <dirent.h> /* readdir */
#include <string.h> /* strlen */
#include <math.h> /* fabs */

/* Serial implementation */
bool invert_matrix(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol]);
bool gaussian_elimination(int nrow, int ncol, double mat[nrow][ncol]);
bool rref(int nrow, int ncol, double mat[nrow][ncol]);

/* Serial helper functions */
void augment_mat(int n, double mat[n][n], double mat_aug[n][2 * n]);
void swap_rows(int r1, int r2, int nrow, int ncol, double mat[nrow][ncol]);
void multiply_row(int row_idx, double n, int nrow, int ncol, double mat[nrow][ncol]);
void subtract_row(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol]);
void extract_inverse(int nrow, int ncol, double mat_aug[nrow][ncol], double mat_inv[nrow][nrow]);

/* Run test matrices */
bool invert_mat_from_file(const char* dir, const char* fname);

/* Testing/debugging functions */
bool check_inverse(int nrow, int ncol, double m1[nrow][ncol], double m2[nrow][ncol]);
void copy_mat(int nrow, int ncol, double src[nrow][ncol], double dest[nrow][ncol]);
void print_mat(int nrow, int ncol, double mat[nrow][ncol]);
void print_working_dir();



int main(int argc, char* argv[]) {
    /* int n = 3; */ /* Size of the matrix */

    /* Get and print input matrix */
    /*
       double mat[3][3] = { {0, 3, 2}, {1, 0, 2}, {0, 0, 1} };
       int nrow = sizeof(mat) / sizeof(mat[0]);
       int ncol = sizeof(mat[0]) / sizeof(mat[0][0]);

       print_mat(nrow, ncol, mat);
       */


    FILE *fp = fopen("HPC.ParallelMatrixInversion/test_matrices/mat_3x3_i_1.txt", "r");
    /* Check that the file exists */
    if (!fp) {
	perror("fopen");
	return 1;
    }	

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

void print_working_dir() {	
    /* Check the working directory */
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
	printf("Current working directory: %s\n", cwd);
    } else {
	perror("getcwd");
    }
}

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

/* A method to run test matrices on invert_matrix */
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
    printf("Reading matrix from file %s\n", fname);
    int nrow, ncol;
    char kind;

    if (sscanf(fname, "mat_%dx%d_%c", &nrow, &ncol, &kind) != 3) {
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
    bool res = invert_matrix(nrow, ncol, mat_cp, mat_inv);

    if (res && kind == 'i') {
	if (check_inverse(nrow, ncol, mat, mat_inv)) {
	    printf("Found correct inverse\n");
	} else {
	    printf("ERROR: Didn't find inverse for invertible matrix\n");
	}
    } else if (!res && kind == 'i') {
	printf("ERROR: Failed to invert invertible matrix\n");
    }

    printf("\n\n");
    return true;
}

bool invert_matrix(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol]) {
    int n = nrow;

    /* Augment identity */
    double mat_aug[n][2 * n];
    augment_mat(n, mat, mat_aug);

    /* Forward elimination */
    bool res = gaussian_elimination(n, 2 * n, mat_aug);
    if (res) {
	printf("GE successsful\n");
    } else {
	printf("GE failed\n");
	return false;
    }

    /*
    printf("Matrix after gaussian elimination (should have non-zero diagonal)\n");
    print_mat(n, 2 * n, mat_aug);
    */

    /* Backward elimination */
    bool res2 = rref(n, 2 * n, mat_aug);
    if (res2) {
	printf("RREF successful\n");
    } else {
	printf("RREF failed\n");
	return false;
    }

    /*
    printf("Matrix after RREF\n");
    print_mat(n, 2 * n, mat_aug);
    */

    /* Extract inverse if it the steps before were successfull */
    extract_inverse(n, 2 * n, mat_aug, mat_inv);

    /*
    printf("The extracted inverse matrix\n");
    print_mat(n, n, mat_inv);
    */

    return true;
}

/* Extract the augmented part of the augmented matrix */
void extract_inverse(int nrow, int ncol, double mat_aug[nrow][ncol], double mat_inv[nrow][nrow]) {
    int i, j;

    for (i = 0; i < nrow; i++) {
	for (j = 0; j < ncol; j++) {
	    /* Copy only the right side of the augmented matrix */
	    mat_inv[i][j] = mat_aug[i][nrow + j];
	}
    }
}

/* Second part of the Gauss-Jordan elimination, results in the reduced row echelon form */
bool rref(int nrow, int ncol, double mat[nrow][ncol]) {
    int i; 
    /* printf("rref input: nrow = %d, ncol = %d, matrix = \n", nrow, ncol); */
    /* print_mat(nrow, ncol, mat); */

    for (i = nrow - 1; i > 0; i--) {
	/* printf("GE: Eliminating the rows above %d\n", i); */
	int r;
	for (r = i-1; r >= 0; r--) {
	    /* printf("GE: Eliminating row %d\n", r); */
	    double coeff = mat[r][i];
	    subtract_row(i, r, coeff, nrow, ncol, mat);
	}
    }
    return true;
}


/* Normalize pivots and clear nonzero values below diagonal */
bool gaussian_elimination(int nrow, int ncol, double mat[nrow][ncol]) {
    int i, r;

    /* Iterate the rows of mat */
    for (i = 0; i < nrow; i++) {
	if (fabs(mat[i][i]) < 1e-9) {


	    /* Find row below the current row where the current column index is nonzero */
	    bool found = false;
	    for (r = i+1; r < nrow; r++) {
		if (fabs(mat[r][i]) > 1e-9) {
		    swap_rows(i, r, nrow, ncol, mat);
		    found = true;
		    break;
		}
	    }
	    if (!found) {
		printf("Matrix is singular or nearly singular\n");
		return false;
	    }

	}

	/* Normalize the row */
	double s = 1/mat[i][i];
	multiply_row(i, s, nrow, ncol, mat);

	/* Eliminate nonzero values below */
	int r;
	for (r = i+1; r < nrow; r++) {
	    double coeff = mat[r][i];
	    subtract_row(i, r, coeff, nrow, ncol, mat);
	}
    }

    return true;
}

/* Subtract the values of row row_idx multiplied with coefficient coeff from the row given by target_idx */
void subtract_row(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol]) {
    /* TODO: debug printing */

    if (row_idx < 0 || row_idx >= nrow) {
	printf("Subtract row: invalid row index %d for matrix %d x %d\n", row_idx, nrow, ncol);
    }

    if (target_idx < 0 || target_idx >= nrow) {
	printf("Subtract row: invalid target row index %d for matrix %d x %d\n", target_idx, nrow, ncol);
    }

    int i;
    for (i = 0; i < ncol; i++) {
	mat[target_idx][i] -= mat[row_idx][i] * coeff;
    }
}

void multiply_row(int row_idx, double s, int nrow, int ncol, double mat[nrow][ncol]) {	
    /* TODO: remove later, random debugging */
    if (row_idx < 0 || row_idx >= nrow) {
	printf("Multiply row: invalid row index %d for matrix %d x %d", row_idx, nrow, ncol);
	return;
    }

    int j;
    for (j = 0; j < ncol; j++) {
	mat[row_idx][j] *= s;
    }
}


void swap_rows(int r1, int r2, int nrow, int ncol, double mat[nrow][ncol]) {
    /* TODO: remove later, random debugging */
    if (r1 < 0 || r1 >= nrow) {
	printf("Swap row: invalid row index %d for %d x %d matrix", r1, nrow, ncol);
    }
    if (r2 < 0 || r2 >= nrow) {
	printf("Swap row: invalid row index %d for %d x %d matrix", r2, nrow, ncol);
    }

    double temp;
    int i;
    for (i = 0; i < ncol; i++) {
	temp = mat[r1][i];
	mat[r1][i] = mat[r2][i];
	mat[r2][i] = temp;
    }
}

/* This surely can be parallelized */
void augment_mat(int n, double mat[n][n], double mat_aug[n][2 * n]) {
    int row, col;

    /* Parallelize the copying of the rows */
    for (row = 0; row < n; row++) {
	for (col = 0; col < n; col++) {
	    /* Copy the row of original matrix */
	    mat_aug[row][col] = mat[row][col];
	    /* Create a row of identity matrix to the right */
	    mat_aug[row][n + col] = (row == col) ? 1 : 0;
	}
    }
}

void print_mat(int nrow, int ncol, double mat[nrow][ncol]) {
    int i, j;

    for (i = 0; i < nrow; i++) {
	for (j = 0; j < ncol; j++) {
	    printf("%.2f \t", mat[i][j]);
	}
	printf("\n");
    }
}

void copy_mat(int nrow, int ncol, double src[nrow][ncol], double dest[nrow][ncol]) {
    int i, j;
    for (i = 0; i < nrow; i++) {
	for (j = 0; j < ncol; j++) {
	    dest[i][j] = src[i][j];
	}
    }
}


