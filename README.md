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

---

## Compiling the Files

## Submitting Jobs to a Cluster

The `matrix_inversion.sh` script is configured with different values for `ncpus` as needed. To submit the script to a cluster, use:

```bash
qsub matrix_inversion.sh
```

This command returns a task ID. The resulting output and errors can be found in:

- `matrix_inversion.sh.o[task_id]`
- `matrix_inversion.sh.e[task_id]`

---

## Folder Structure

- **Source Files**:
  - `main.c`: Main file.
  - `common`: Utility file.
- **Shell Script**: `matrix_inversion.sh` for submitting cluster jobs.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---
