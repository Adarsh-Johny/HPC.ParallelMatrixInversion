#!/bin/bash
# Configuration for OpenMP
#PBS -l select=1:ncpus=8:mem=2gb
#PBS -l walltime=0:00:10
#PBS -q short_cpuQ

# load the module
module load openmpi-4.0.4

# Actual running command
./HPC.ParallelMatrixInversion/matrix_inversion

