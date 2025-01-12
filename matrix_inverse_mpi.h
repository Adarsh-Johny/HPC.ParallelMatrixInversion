#ifndef MPI_MATRIX_INVERSE_H
#define MPI_MATRIX_INVERSE_H

void inverse_matrix_mpi(double **mat, int nrow, int ncol, double **mat_inv_parallel);
void benchmark_inversion(double **mat, int nrow, int ncol);

#endif // MPI_MATRIX_INVERSE_H
