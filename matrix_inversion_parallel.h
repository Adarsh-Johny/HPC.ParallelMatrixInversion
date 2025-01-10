#ifndef MATRIX_INVERSION_PARALLEL_H
#define MATRIX_INVERSION_PARALLEL_H
#include <stdbool.h>

bool invert_matrix_par(int nrow, int ncol, double mat[nrow][ncol], double mat_inv[nrow][ncol]);
bool gaussian_elimination_par(int nrow, int ncol, double mat[nrow][ncol]);
bool rref_par(int nrow, int ncol, double mat[nrow][ncol]);
void extract_inverse_par(int nrow, int ncol, double mat_aug[nrow][ncol], double mat_inv[nrow][nrow]);
void augment_mat_par(int n, double mat[n][n], double mat_aug[n][2 * n]);
void subtract_row_par(int row_idx, int target_idx, double coeff, int nrow, int ncol, double mat[nrow][ncol]);
void multiply_row_par(int row_idx, double s, int nrow, int ncol, double mat[nrow][ncol]);

void benchmark_matrix_inversion(int nrow, int ncol, double mat[nrow][ncol]);

#endif
