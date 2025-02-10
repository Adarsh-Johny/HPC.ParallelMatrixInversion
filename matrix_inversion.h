#ifndef MATRIX_INVERSION_H
#define MATRIX_INVERSION_H

#include <stdbool.h>

/* Implementation */
double** invert_matrix(int nrow, int ncol, double** mat, bool *success);
bool gaussian_elimination(int nrow, int ncol, double** mat);
bool rref(int nrow, int ncol, double** mat);

/* Helper functions */
void augment_mat(int n, double** mat, double** mat_aug);
void multiply_row(int row_idx, double n, int nrow, int ncol, double** mat);
void subtract_row(int row_idx, int target_idx, double coeff, int nrow, int ncol, double** mat);
double** extract_inverse(int nrow, int ncol, double** mat_aug);

#endif
