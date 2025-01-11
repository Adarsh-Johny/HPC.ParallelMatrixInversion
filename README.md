# HPC.ParallelMatrixInversion
Parallel Matrix Inversion using MPI for High Performance Computing - Data Science



gcc -std=c99 -g -Wall -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c main.c -lm
