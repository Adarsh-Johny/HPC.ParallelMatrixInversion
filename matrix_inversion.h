#ifndef MATRIX_INVERSION_H
#define MATRIX_INVERSION_H

#include <stdbool.h>

/* Serial implementation */
bool invert_matrix(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol]);
bool gaussian_elimination(int nrow, int ncol, double mat[nrow][ncol]);
bool rref(int nrow, int ncol, double mat[nrow][ncol]);

/* Serial helper functions */
void augment_mat(int n, double mat[n][n], double mat_aug[n][2 * n]);
void multiply_row(int row_idx, double n, int nrow, int ncol, double mat[nrow][ncol]);
void subtract_row(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol]);
void extract_inverse(int nrow, int ncol, double mat_aug[nrow][ncol], double mat_inv[nrow][nrow]);

/* Function for benchmarking */
void benchmark_matrix_inversion(int nrow, int ncol, double mat[nrow][ncol]);

#endif