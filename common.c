#include <stdio.h> /* NULL */
#include <unistd.h> /* getcwd */

/* Swap rows r1 and r2 */
void swap_rows(int r1, int r2, int nrow, int ncol, double mat[nrow][ncol]) {
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
