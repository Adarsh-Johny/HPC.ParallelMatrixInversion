# HPC.ParallelMatrixInversion

In this project, we implemented a parallel matrix inversion algorithm using OpenMP and MPI as part of the High-Performance Computing course at the University of Trento.

---

## Prerequisites

To run this project, ensure the following packages and tools are installed:

### System Packages

- **C Compiler**: `gcc` or `mpicc` with support for C99 or later.
- **MPI Library**: `OpenMPI`.
- **Python**: Python 3.x for Jupyter Notebook.

#### Install Necessary Packages

**Ubuntu/Debian:**

```bash
sudo apt update
sudo apt install build-essential libopenmpi-dev python3 python3-pip
```

**Fedora:**

```bash
sudo dnf install gcc openmpi openmpi-devel python3 python3-pip
```

**Arch Linux:**

```bash
sudo pacman -S gcc openmpi python python-pip
```

### Python Packages for Jupyter Notebook

The Jupyter Notebook (`matrix_generator.ipynb`) requires the following Python packages:

- `numpy`
- `matplotlib`
- `jupyter`

Install them using pip:

```bash
pip install numpy matplotlib jupyter
```

---

## Description of Implementations

This project includes three implementations of the matrix inversion algorithm:

1. **OpenMP Implementation**
   - Main File: `main.c`
   - Utilizes OpenMP for parallel computation.

2. **MPI Implementation**
   - Main File: `mpi_inverse_main.c`
   - Leverages MPI for distributed-memory parallelism.

3. **Serial Implementation**
   - Main File: `main_serial.c`
   - Provides a baseline for performance comparison.

---

## Compiling the Files

1. **OpenMP Execution** (Main File: `main.c`)

   ```bash
   mpicc -std=c99 -g -Wall -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c main.c -lm
   ```

2. **MPI Execution** (Main File: `mpi_inverse_main.c`)

   ```bash
   mpicc -std=c99 -g -Wall -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c mpi_inverse_main.c -lm
   ```

3. **Serial Execution** (Main File: `main_serial.c`)

   ```bash
   mpicc -std=c99 -g -Wall -fopenmp -I./helpers -o main_program ./helpers/common.c ./helpers/file_reader.c matrix_inversion_parallel.c matrix_inversion.c main_serial.c -lm -pg
   ```

---

## Submitting Jobs to a Cluster

The `matrix_inversion.sh` script is configured with different values for `ncpus` as needed. To submit the script to a cluster, use:

```bash
qsub matrix_inversion.sh
```

This command returns a task ID. The resulting output and errors can be found in:

- `matrix_inversion.sh.o[task_id]`
- `matrix_inversion.sh.e[task_id]`

---

## Generating Test Matrices and Performance Metrics

- **Copy paste the logs and rename them based on the naming specified in the jupiter notebook**
- **Matrix Generator**: Use the `matrix_generator.ipynb` notebook to generate test matrices.
- **Performance Metrics**: Metrics are automatically saved in the `Metrics` folder after code execution.

---

## Folder Structure

- **Source Files**:
  - `main.c`: OpenMP implementation.
  - `mpi_inverse_main.c`: MPI implementation.
  - `main_serial.c`: Serial implementation.
  - `helpers/`: Contains utility files (`common.c`, `file_reader.c`).
- **Jupyter Notebook**: `matrix_generator.ipynb` for generating test matrices.
- **Metrics Folder**: Stores performance metrics.
- **Shell Script**: `matrix_inversion.sh` for submitting cluster jobs.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---
