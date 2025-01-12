#ifndef MATRIX_INVERSION_H
#define MATRIX_INVERSION_H

#include <stdbool.h>

/* Serial implementation */
bool invert_matrix(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol]);
bool gaussian_elimination(int nrow, int ncol, double mat[nrow][ncol]);
bool rref(int nrow, int ncol, double mat[nrow][ncol]);

/* Serial helper functions */
void augment_mat_ser(int n, const double mat[n][n], double mat_aug[n][2 * n]);
void extract_inverse_ser(int nrow, int ncol, const double mat_aug[nrow][ncol], double mat_inv[nrow][nrow]);
void multiply_row_ser(int row_idx, double n, int nrow, int ncol, double mat[nrow][ncol]);
void subtract_row_ser(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol]);

/* Function for benchmarking */
void benchmark_matrix_inversion(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol]);

#endif