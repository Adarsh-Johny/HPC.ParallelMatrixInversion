#!/bin/bash
# #PBS -N generate_matrix
#PBS -l select=1:ncpus=4:mem=5GB
#PBS -l walltime=00:10:00
#PBS -q short_cpuQ
# #PBS -o output.log
# #PBS -e error.log

source /etc/profile.d/modules.sh  

# Load a Python module and hope it contains numpy
module load python-3.7.2

# Check if it does
python3 -c "import numpy; print('numpy is available')"

module load python-3.7.2
python3 ./HPC.ParallelMatrixInversion/generate_matrix.py

