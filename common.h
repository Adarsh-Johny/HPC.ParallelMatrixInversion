#ifndef COMMON_H
#define COMMON_H

void copy_mat(int nrow, int ncol, double** src, double** dest);
void print_mat(int nrow, int ncol, double** mat);
void print_working_dir();
void swap_rows(int r1, int r2, int nrow, int ncol, double** mat);
double** allocate_mat(int nrow, int ncol);
void free_mat(double** mat, int nrow);

#endif
