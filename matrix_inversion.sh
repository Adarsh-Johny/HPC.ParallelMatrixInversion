#!/bin/bash
#PBS -l select=1:ncpus=16:mem=2gb
#PBS -l walltime=0:05:00
#PBS -q short_cpuQ

# Load the module
module load openmpi-4.0.4

# Set the number of OpenMP threads
# export OMP_NUM_THREADS=4

export OMP_DISPLAY_ENV=TRUE

# Define file paths
MATRIX_FILES=(
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_5x5_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_8x8_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_10x10_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_20x20_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_30x30_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_40x40_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_50x50_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_60x60_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_70x70_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_80x80_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_90x90_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_100x100_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_200x200_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_300x300_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_400x400_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_500x500_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_600x600_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_700x700_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_1000x1000_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_1500x1500_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_2000x2000_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_2500x2500_01.txt"
              "HPC.ParallelMatrixInversion/performance_test_matrices/matrix_3000x3000_01.txt"
              )

# Loop through files and execute the program for each
for FILE_PATH in "${MATRIX_FILES[@]}"; do
    echo "Processing file: $FILE_PATH"
    mpiexec ./HPC.ParallelMatrixInversion/main_program -path="$FILE_PATH"

done