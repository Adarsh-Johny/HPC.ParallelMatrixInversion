#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

void print_mat(int nrow, int ncol, double mat[nrow][ncol]);

int main(int argc, char* argv[]) {
	
	double mat[2][2] = { {1, 2}, {0, 1} };
	int nrow = sizeof(mat) / sizeof(mat[0]);
	int ncol = sizeof(mat[0]) / sizeof(mat[0][0]);

	print_mat(nrow, ncol, mat);

	return 0;
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
