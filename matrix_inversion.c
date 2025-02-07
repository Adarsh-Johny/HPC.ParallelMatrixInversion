#include "matrix_inversion.h"
#include "common.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h> /* fabs */

bool invert_matrix(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol]) {
    int n = nrow;

    /* Augment identity */
    double mat_aug[n][2 * ncol];
    augment_mat(n, mat, mat_aug);

    /* Forward elimination */
    bool res = gaussian_elimination(n, 2 * n, mat_aug);
    if (res) {
	//printf("GE successsful\n");
    } else {
	//printf("GE failed\n");
	return false;
    }

    /*
    printf("Matrix after gaussian elimination (should have non-zero diagonal)\n");
    print_mat(n, 2 * n, mat_aug);
    */

    /* Backward elimination */
    bool res2 = rref(n, 2 * n, mat_aug);
    if (res2) {
	//printf("RREF successful\n");
    } else {
	printf("ERROR: RREF failed\n");
	return false;
    }

    /*
    printf("Matrix after RREF\n");
    print_mat(n, 2 * n, mat_aug);
    */

    /* Extract inverse */
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


