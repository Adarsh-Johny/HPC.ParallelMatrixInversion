#include <stdio.h> /* NULL */
#include <stdlib.h> /* free */
#include <unistd.h> /* getcwd */

void free_mat(double** mat, int nrow) {
    if (mat == NULL) {
	printf("Free_mat: tried to free a NULL pointer\n");
	return;
    }

    int i;
    for (i = 0; i < nrow; i++) {
	free(mat[i]);
    }
    free(mat);
}

double** allocate_mat(int nrow, int ncol) {
    double** mat = malloc(nrow * sizeof(double));
    int i;
    for (i = 0; i < nrow; i++) {
	mat[i] = malloc(ncol * sizeof(double));
    }
    return mat;
}

/* Swap rows r1 and r2 */
void swap_rows(int r1, int r2, int nrow, int ncol, double** mat) {
    double temp;
    int i;
    for (i = 0; i < ncol; i++) {
	temp = mat[r1][i];
	mat[r1][i] = mat[r2][i];
	mat[r2][i] = temp;
    }
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

void print_mat(int nrow, int ncol, double** mat) {
    int i, j;

    for (i = 0; i < nrow; i++) {
	for (j = 0; j < ncol; j++) {
	    printf("%.2f \t", mat[i][j]);
	}
	printf("\n");
    }
}

void copy_mat(int nrow, int ncol, double** src, double** dest) {
    int i, j;
    for (i = 0; i < nrow; i++) {
	for (j = 0; j < ncol; j++) {
	    dest[i][j] = src[i][j];
	}
    }
}
