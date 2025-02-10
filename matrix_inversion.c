#include "matrix_inversion.h"
#include "common.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h> /* fabs */

double** invert_matrix(int nrow, int ncol, double** mat, bool *success) {
    int n = nrow;
    *success = false;

    /* Augment identity */
    double** mat_aug = allocate_mat(nrow, 2 * ncol);
    if (!mat_aug) {
	printf("Error in invert_matrix: memory allocation failed\n");
	return NULL;
    }

    augment_mat(n, mat, mat_aug);

    /* Forward elimination */
    if (!gaussian_elimination(n, 2 * n, mat_aug)) {
	//printf("GE failed\n");
	free_mat(mat_aug, n);
	return NULL;
    }

    /*
    printf("Matrix after gaussian elimination (should have non-zero diagonal)\n");
    print_mat(n, 2 * n, mat_aug);
    */

    /* Backward elimination */
    if (!rref(n, 2 * n, mat_aug)) {
	printf("ERROR: RREF failed\n");
	free_mat(mat_aug, n);
	return NULL;
    }

    /*
    printf("Matrix after RREF\n");
    print_mat(n, 2 * n, mat_aug);
    */

    /* Extract inverse */
    double** mat_inv = extract_inverse(n, 2 * n, mat_aug);

    /*
    printf("The extracted inverse matrix\n");
    print_mat(n, n, mat_inv);
    */

    free_mat(mat_aug, n);
    *success = true;
    return mat_inv;
}

/* Extract the augmented part of the augmented matrix */
double** extract_inverse(int nrow, int ncol, double** mat_aug) {
    double** mat_inv = allocate_mat(nrow, nrow);
    if (!mat_inv) {
	printf("Error in extract_inverse: memory allocation failed\n");
	return NULL;
    }

    int i, j;

    //#pragma omp parallel for
    for (i = 0; i < nrow; i++) {
	for (j = 0; j < nrow; j++) {
	    /* Copy only the right side of the augmented matrix */
	    mat_inv[i][j] = mat_aug[i][nrow + j];
	}
    }
    return mat_inv;
}

/* Second part of the Gauss-Jordan elimination, results in the reduced row echelon form */
bool rref(int nrow, int ncol, double** mat) {
    int i; 

    /* Eliminate the values above the pivots starting from the last row */
    for (i = nrow - 1; i > 0; i--) {
	int row;
        //#pragma omp parallel for // No because of the coeff somehow gets messed up
	for (row = i - 1; row >= 0; row--) {
	    double coeff = mat[row][i];
	    int col;
	    for (col = 0; col < ncol; col++) {
		mat[row][col] -= mat[i][col] * coeff;
	    }
	}
    }
    return true;
}


/* Normalize pivots and clear nonzero values below diagonal */
bool gaussian_elimination(int nrow, int ncol, double** mat) {
    int i;

    /* Iterate the rows of mat */
    for (i = 0; i < nrow; i++) {
	if (fabs(mat[i][i]) < 1e-9) {


	    /* Find row below the current row where the current column index is nonzero */
	    bool found = false;
	    int row;
	    for (row = i + 1; row < nrow; row++) {
		if (fabs(mat[row][i]) > 1e-9) {
		    swap_rows(i, row, nrow, ncol, mat);
		    found = true;
		    break;
		}
	    }
	    if (!found) {
		printf("Matrix is singular or nearly singular\n");
		return false;
	    }

	}

	/* Normalize the pivot to 1 by scaling the row */
	double s = 1 / mat[i][i];
	int col;
	for (col = 0; col < ncol; col++) {
	    mat[i][col] *= s;
	}

	/* Eliminate nonzero values below the current pivot */
	int row;
	for (row = i + 1; row < nrow; row++) {
	    double coeff = mat[row][i];
	    int col;
	    for (col = 0; col < ncol; col++) {
		mat[row][col] -= coeff * mat[i][col];
	    }
	}
    }

    return true;
}

void augment_mat(int n, double** mat, double** mat_aug) {
    int row, col;

    //#pragma omp parallel for
    for (row = 0; row < n; row++) {
        //#pragma omp parallel for
	for (col = 0; col < n; col++) {
	    /* Copy the row of original matrix */
	    mat_aug[row][col] = mat[row][col];
	    /* Create a row of identity matrix to the right */
	    mat_aug[row][n + col] = (row == col) ? 1 : 0;
	}
    }
}


