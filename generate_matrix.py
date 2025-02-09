import numpy as np

def generate_invertible_matrix(n, seed=None):
    if seed is not None:
        np.random.seed(seed)

    # Keeping the value range small might help with numerical problems for bigger matrices
    min_val = -2
    max_val = 2

    while True:
        A = np.random.rand(n, n) * (max_val-min_val) + min_val

        if np.linalg.cond(A) < 1 / np.finfo(A.dtype).eps: # check the condition number
            print(np.finfo(A.dtype).eps)
        #if np.linalg.det(A) != 0:
            return A

mat_sizes = [2**5, 2**6, 2**7, 2**8, 2**9, 2**10, 2**11, 2**12,]# 2**13]

for size in mat_sizes:
    mat = generate_invertible_matrix(size)
    np.savetxt(f"./HPC.ParallelMatrixInversion/data/mat_{size}.txt", mat)
    print(f"Saved matrix of size {size}")

