#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

void print_mat(int nrow, int ncol, double mat[nrow][ncol]);
void augment_mat(int n, double mat[n][n], double mat_aug[n][2 * n]);
void swap_rows(int r1, int r2, int nrow, int ncol, double mat[nrow][ncol]);
void multiply_row(int row_idx, double n, int nrow, int ncol, double mat[nrow][ncol]);
void subtract_row(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol]);
int gaussian_elimination(int nrow, int ncol, double mat[nrow][ncol]);
int rref(int nrow, int ncol, double mat[nrow][ncol]);
void extract_inverse(int nrow, int ncol, double mat_aug[nrow][ncol], double mat_inv[nrow][nrow]);
int invert_matrix(int nrow, int ncol, double mat[nrow][ncol]);

int main(int argc, char* argv[]) {
	/* int n = 3; */ /* Size of the matrix */
	
	/* Get and print input matrix */
	double mat[3][3] = { {0, 3, 2}, {1, 0, 2}, {0, 0, 1} };
	int nrow = sizeof(mat) / sizeof(mat[0]);
	int ncol = sizeof(mat[0]) / sizeof(mat[0][0]);

	print_mat(nrow, ncol, mat);

	int success = invert_matrix(nrow, ncol, mat);
	printf("Inversion success: %d", success);

	/* If not, augment matrix and print the result */
	/* 
	* double mat_aug[n][2 * n];
	* augment_mat(n, mat, mat_aug);
	*
	* printf("Augmented matrix:\n");
	* print_mat(n, 2 * n, mat_aug);
	*/

	/*
	* swap_rows(0, 1, n, 2 * n, mat_aug);
	* 
	* printf("Augmented matrix after row swap\n");
	* print_mat(n, 2 * n, mat_aug);
	* */

	/* multiply_row(1, 5, n, 2 * n, mat_aug); */

	/* printf("Augmented matrix after row multiplication\n"); */
	/* print_mat(n, 2 * n, mat_aug); */
	/*
	* subtract_row(1, 2, 1, n, 2 * n, mat_aug);
	*
	* printf("Augmented matrix after row subtraction\n");
	* print_mat(n, 2 * n, mat_aug);
	* */

	return 0;
}

int invert_matrix(int nrow, int ncol, double mat[nrow][ncol]) {
	int n = nrow;

	/* Augment identity */
	double mat_aug[n][2 * n];
	augment_mat(n, mat, mat_aug);

	/* Forward elimination */
	int res = gaussian_elimination(n, 2 * n, mat_aug);
	if (res == 1) {
		printf("GE successsful\n");
	} else {
		printf("GE failed\n");
		return 0;
	}

	printf("Matrix after gaussian elimination (should have non-zero diagonal)\n");
	print_mat(n, 2 * n, mat_aug);

	
	/* Backward elimination */
	int res2 = rref(n, 2 * n, mat_aug);
	if (res2 == 1) {
		printf("RREF successful\n");
	} else {
		printf("RREF failed\n");
		return 0;
	}

	printf("Matrix after RREF\n");
	print_mat(n, 2 * n, mat_aug);

	/* Extract inverse */
	double mat_inv[n][n];
	extract_inverse(n, 2 * n, mat_aug, mat_inv);

	printf("The extracted inverse matrix\n");
	print_mat(n, n, mat_inv);
	
	return 1;
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
int rref(int nrow, int ncol, double mat[nrow][ncol]) {
	int i; 
	printf("rref input: nrow = %d, ncol = %d, matrix = \n", nrow, ncol);
	print_mat(nrow, ncol, mat);

	for (i = nrow - 1; i > 0; i--) {
		printf("GE: Eliminating the rows above %d\n", i);
		int r;
		for (r = i-1; r >= 0; r--) {
			printf("GE: Eliminating row %d\n", r);
			double coeff = mat[r][i];
			subtract_row(i, r, coeff, nrow, ncol, mat);
		}
	}
	return 1;
}


/* Results in unreduced row echelon form in which the pivots are normalized.
 * Operates on a square matrix or the square part of an augmented matrix. */
int gaussian_elimination(int nrow, int ncol, double mat[nrow][ncol]) {
	int i, r;

	/* Iterate the rows of mat */
	for (i = 0; i < nrow; i++) {
		printf("GE: i = %d, diagonal value %f\n", i, mat[i][i]);
		
		/* Check if the diagonal is zero: if it is, swap find another row where the current column 
		 * contains a non-zero value */
		if (mat[i][i] == 0) {
			/* If we are on the last row and the diagonal is 0, 
			 * there's a zero column. The matrix is not invertible so return 0 */
			if (i == nrow - 1) {
				printf("GE: there's something wrong with the matrix\n");
				print_mat(nrow, ncol, mat);
				return 0;
			}

			/* Find row below the current row in which there's no 0 */
			for (r = i+1; r < nrow; r++) {
				if (mat[r][i] != 0) {
					swap_rows(i, r, nrow, ncol, mat);
					printf("GE: rows %d and %d swapped, resulting matrix\n", i, r);
					print_mat(nrow, ncol, mat);
				}
			}

		} else {
			/* Normalize and eliminate */
			double s = 1/mat[i][i];
			multiply_row(i, s, nrow, ncol, mat);
			printf("Row %d normalized\n", i);

			/* If we are on the last row, finish */
			if (i == nrow - 1) {
				printf("GE: We are on the last row %d/%d\n", i, nrow);
				return 1;
			}

			printf("GE: Eliminating the rows below %d\n", i);
			int r;
			for (r = i+1; r < nrow; r++) {
				printf("GE: Eliminating row %d\n", r);
				double coeff = mat[r][i];
				subtract_row(i, r, coeff, nrow, ncol, mat);				
			}

			return 1;
		}
	}
	printf("GE: diagonal is nonzero\n");
	return 1;
}

/* Subtract the values of row row_idx multiplied with coefficient coeff from the row given by target_idx */
void subtract_row(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol]) {
	/* TODO: debug printing */
	printf("Subtracting %d * %.2f from %d\n", row_idx, coeff, target_idx);

	if (row_idx < 0 || row_idx >= nrow) {
		printf("Subtract row: invalid row index %d for matrix %d x %d\n", row_idx, nrow, ncol);
	}
	
	if (target_idx < 0 || target_idx >= nrow) {
		printf("Subtract row: invalid target row index %d for matrix %d x %d\n", target_idx, nrow, ncol);
	}

	int i;
	for (i = 0; i < ncol; i++) {
		printf("i = %d\n", i);
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


/* Possibly parallelizable or to combine with sth else */
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


