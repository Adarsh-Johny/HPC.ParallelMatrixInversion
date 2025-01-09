#ifndef COMMON_H
#define COMMON_H

void copy_mat(int nrow, int ncol, double src[nrow][ncol], double dest[nrow][ncol]);
void print_mat(int nrow, int ncol, double mat[nrow][ncol]);
void print_working_dir();
void swap_rows(int r1, int r2, int nrow, int ncol, double mat[nrow][ncol]);

#endif
