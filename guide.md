## Connect to the cluster

- Turn GlobalProtect VPN on
- `ssh susanna.repo@hpc2.unitn.it` (or `ssh -vvvv susanna.repo@hpc2.unitn.it` to show debugging)
- `vpn-out.icts.unitn.it` - VPN URL

## Loading modules

- List the loaded modules: `module list`
- List the available modules: `module avail`

# Using mpich-3.2 library

## Load the mpich module

- List the loaded modules: `module list` and see if `mpich-3.2` is loaded
  - If not, `module load mpich-3.2`
- List the available modules: `module avail` if something else is missing

## Copy a file from local machine to the cluster

Run `scp test_code.c susanna.repo@hpc2.unitn.it:/home/susanna.repo`

## Write a script to run `mpi_code` using mpich-3.2 on the cluster

- Create the shell script, e.g. `test_script.sh`:

```sh
#!/bin/bash
# 2 nodes, 10 cpus/node, 2 GB RAM
#PBS -l select=2:ncpus=10:mem=2gb
# set max exec time (1 min)
#PBS -l walltime=0:01:00
# Put the job in the short queue (max. 6 h)
#PBS -q short_cpuQ
# Set MPI environment
module load mpich-3.2
# Actual running command using 20 processes
mpirun.actual -n 20 ./mpi_code
```

## Running MPI programs with mpich

- Write the source: `vi mpi_code.c`
- Compile the source code: `mpicc -g -Wall -o mpi_code mpi_code.c`
  - `-g` for debugging
  - `-Wall` to show warnings
  - `-o` to default the output file name `mpi_code`
  - Add `-std=c11` to use the newer C standard if something causes an error
- Edit the script: set nodes, cpus and memory etc.
- Submit the job: `qsub test_script.sh`
  - You get back the job ID!
- Check the status of the job: `qstat -f <job_id>`, (`-f` not needed) more info with `tracejob <job_id>`
  - Cancel the job: `qdel <job_id>`
- Check the new err and out files: `test_script.sh.e` and `test_script.sh.o`:
  - e.g. `more test_script.sh.e`, prints the errors to the console

# OpenMPI

## Load modules

- `module load openmpi-4.0.4`

## Running MPI programs with openMPI in the cluster

- Write the source
- Compile it by running `gcc −g −Wall −fopenmp −o omp_hello omp_hello.c`
  - Without cluster: run with 4 threads (passed as an argc to main): `./omp_hello 4`
- Create script file, for example (without args):

    ```sh
    #!/bin/bash
    #PBS -l select=2:ncpus=10:mem=2gb
    #PBS -l walltime=0:01:00
    #PBS -q short_cpuQ

    # load the module
    module load openmpi-4.0.4

    # Actual running command
    ./exercises-openmp/hello_world
    ```
