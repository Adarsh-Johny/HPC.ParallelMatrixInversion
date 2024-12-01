#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

void print_mat(int nrow, int ncol, double mat[nrow][ncol]);
void augment_mat(int n, double mat[n][n], double aug_mat[n][2 * n]);
void swap_rows(int r1, int r2, int nrow, int ncol, double mat[nrow][ncol]);
void multiply_row(int row_idx, double n, int nrow, int ncol, double mat[nrow][ncol]);

int main(int argc, char* argv[]) {
	
	/* Get and print input matrix */
	double mat[2][2] = { {1, 2}, {0, 1} };
	int nrow = sizeof(mat) / sizeof(mat[0]);
	int ncol = sizeof(mat[0]) / sizeof(mat[0][0]);

	print_mat(nrow, ncol, mat);
	
	/* Check if the matrix contains 0-rows or 0-columns */
	

	/* If not, augment matrix and print the result */
	double aug_mat[2][2 * 2];

	augment_mat(2, mat, aug_mat);
	
	printf("Augmented matrix:\n");
	print_mat(2, 2 * 2, aug_mat);

	swap_rows(0, 1, 2, 2 * 2, aug_mat);

	printf("Augmented matrix after row swap\n");
	print_mat(2, 2 * 2, aug_mat);

	multiply_row(1, 5, 2, 2 * 2, aug_mat);

	printf("Augmented matrix after row multiplication\n");
	print_mat(2, 2 * 2, aug_mat);

	return 0;
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
	double temp;
	int i;
	for (i = 0; i < ncol; i++) {
		temp = mat[r1][i];
		mat[r1][i] = mat[r2][i];
		mat[r2][i] = temp;
	}
}

/* This surely can be parallelized */
void augment_mat(int n, double mat[n][n], double aug_mat[n][2 * n]) {
	int row, col;

	/* Parallelize the copying of the rows */
	for (row = 0; row < n; row++) {
		for (col = 0; col < n; col++) {
			/* Copy the row of original matrix */
			aug_mat[row][col] = mat[row][col];
			/* Create a row of identity matrix to the right */
			aug_mat[row][col + n] = (row == col) ? 1 : 0;
		}
	}
}

void print_mat(int nrow, int ncol, double mat[nrow][ncol]) {
	int i, j;

	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			printf("%f \t", mat[i][j]);
		}
		printf("\n");
	}
}
