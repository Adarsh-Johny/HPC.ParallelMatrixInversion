# HPC.ParallelMatrixInversion

In this project, we implemented a parallel matrix inversion algorithm using OpenMPI and MPI as a part of High Performance Computing course in University of Trento.

## Compiling the files

1. OpenMP execution (main file: `main.c`)

`mpicc -std=c99 -g -Wall -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c main.c -lm`

2. MPI (main file: `mpi_inverse_main.c`)

`mpicc -std=c99 -g -Wall -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c mpi_inverse_main.c -lm`

3. Serial (main file: `main_serial.c`)

`mpicc -std=c99 -g -Wall -fopenmp -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c main_serial.c -lm -pg`

`matrix_inversion.sh` is modified according to need with ncpus of different values.

This bash file is then submitted to the cluster by `qsub matrix_inversion.sh` which returns the task id. The resulting output and errors can be read from `matrix_inversion.sh.o[task_id]` and `matrix_inversion.sh.e[task_id]`.

Running the code will also generate performance metrics in the `Metrics` folder. For this, we use the generated test matrices which is done in `matrix_generator.ipynb`.
