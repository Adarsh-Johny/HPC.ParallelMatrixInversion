#ifndef COMMON_H
#define COMMON_H

void copy_mat(int nrow, int ncol, double src[nrow][ncol], double dest[nrow][ncol]);

void print_mat(int nrow, int ncol, double mat[nrow][ncol]);

void print_working_dir();

void swap_rows(int r1, int r2, int nrow, int ncol, double mat[nrow][ncol]);

bool check_inverse(int nrow, int ncol, double m1[nrow][ncol], double m2[nrow][ncol]);

void copy_matrix(int nrow, int ncol, double **source, double dest[nrow][ncol]);

bool compare_matrices(int nrow, int ncol, double mat1[nrow][ncol], double mat2[nrow][ncol]);

void compare_inversions(const char *fname, int nrow, int ncol, double mat_inv_serial[nrow][ncol], double mat_inv_parallel[nrow][ncol]);

void free_matrix(double **mat, int nrow);

#endif