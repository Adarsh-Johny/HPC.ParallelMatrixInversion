# HPC.ParallelMatrixInversion

Parallel Matrix Inversion using MPI for High Performance Computing - Data Science

Compile files:

1. OpenMP execution(main file: `main.c`)

`mpicc -std=c99 -g -Wall -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c main.c -lm`

2. MPI(main file: `mpi_inverse_main.c`)

`mpicc -std=c99 -g -Wall -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c mpi_inverse_main.c -lm`

3. Serial(main file: `main_serial.c`)

`mpicc -std=c99 -g -Wall -fopenmp -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c main_serial.c -lm -pg`

`matrix_inversion.sh` is modified according to need with ncpus of different values

Execute in cluster by this command:
`qsub matrix_inversion.sh`
