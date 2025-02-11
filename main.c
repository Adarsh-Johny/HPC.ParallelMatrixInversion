#include "matrix_inversion.h"
#include "matrix_inversion_parallel.h"
#include "common.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h> /* malloc */
#include <math.h> /* fabs */
#include <dirent.h> /* DIR, readdir */
#include <string.h> /* strlen */


/* Run test matrices */
bool invert_mat_from_file(const char* dir, const char* fname);
bool check_inverse(int nrow, int ncol, double** m1, double** m2);
int run_test_matrices();
void save_results(const char* fname, int size, int threads, int run, double time);
bool run_benchmark_for_file(const char* dir, const char* fname, int repeats, int num_threads, const char* save_location);
void run_benchmark(const char* save_location, int repeats);

int main(int argc, char* argv[]) {

    run_benchmark("./HPC.ParallelMatrixInversion/results_serial.csv", 5);

    /*
    int repeats = 5;
    int i;
    for (i = 0; i < repeats; i++) {
	printf("----------------------");
	printf("Repeat %d", i);
	printf("----------------------\n");
        //int failed = run_test_matrices();
	//printf("Failed tests: %d\n", failed);
	//run_benchmark("./HPC.ParallelMatrixInversion/results_8192.csv");
	printf("----------------------");
	printf("-------");
	printf("----------------------\n\n");
    }
    */
    return 0;
}

void run_benchmark(const char* save_location, int repeats) {
    //printf("DEBUG: run_benchmark\n");
    //int fails = run_test_matrices;

    //if (fails > 0) {
    //    printf("%d tests failed");
    //}
    
    int num_threads;
    #ifdef _OPENMP
        num_threads = omp_get_max_threads();
    #else
        num_threads = 1;
    #endif

    printf("Running with %d threads\n", num_threads);

    print_working_dir();

    const char *directory_path = "./HPC.ParallelMatrixInversion/data/";
    DIR *dir = opendir(directory_path);
    if (!dir) {
            perror("opendir");
	    return;
    }

    printf("Opened dir %s\n", directory_path);

    // Create the file to save the results (we want to create new result file every time)
    FILE *results_file = fopen(save_location, "w");
    if (!results_file) {
	    perror("Error creating results file");
	        return;
    }
    fclose(results_file);

    //printf("DEBUG: Created the file for results\n");
    /* Iterate the matrix data in data folder and time inversion of each matrix as many
     * times as indicated by repeats */
    //int repeats = 2;
    struct dirent *entry;
    //printf("DEBUG: reading from dir %s\n", directory_path);
    while ((entry = readdir(dir)) != NULL) {
	//printf("DEBUG: Found file %s\n", entry->d_name);
        /* Skip directories (there shouldn't be any) */
	//printf("DEBUG: entry type is %d, expected %d\n", entry->d_type, DT_REG);
        if (entry->d_type == DT_REG) {
            printf("Running benchmark for file %s", entry->d_name);
            run_benchmark_for_file(directory_path, entry->d_name, repeats, num_threads, save_location);
	}
    }
    closedir(dir);
}

/* Run benchmark for a given file.
 * The matrix is read from the file, then it's inverted with timing.
 * Inversion is repeated as many times as given in the argument repeats.
 * Each result is saved into same file for each matrix, containing the
 * time, matrix size, number of threads and index.
 */
bool run_benchmark_for_file(const char* dir, const char* fname, int repeats, int num_threads, const char* save_location) {

    char* full_path = malloc(strlen(dir) + strlen(fname) + 1);
    if (!full_path) {
        perror("malloc (concatenate file path)");
        return false;
    }

    strcpy(full_path, dir);
    strcat(full_path, fname);

    int n;

    if (sscanf(fname, "mat_%d", &n) != 1) {
        printf("Invalid filename format: %s\n", fname);
	free(full_path);
        return false;
    }

    FILE *fp = fopen(full_path, "r");
    if (!fp) {
        perror("Error opening file");
	free(full_path);
        print_working_dir();
        return false;
    }

    int i, j;
    //double mat[n][n];
    double** mat = allocate_mat(n, n);
    if (!mat) {
	printf("Error in run_benchmark_for_file: memory allocation (mat) failed");
	free(full_path);
	return false;
    }

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (fscanf(fp, "%lf", &mat[i][j]) != 1) {
                printf("Error reading matrix value at [%d][%d] in file: %s\n", i, j, fname);
                fclose(fp);
		//free_mat(mat, n);
		free(full_path);
                return false;
            }
        }
    }

    fclose(fp);
    free(full_path);
    
    //double mat_cp[n][n];
    double** mat_cp = allocate_mat(n, n);
    if (!mat_cp) {
	printf("Error in run_benchmark_for_file: matrix allocation (mat_cp)  failed\n");
	free_mat(mat, n);
	return false;
    }

    copy_mat(n, n, mat, mat_cp);

    //printf("DEBUG: run_benchmark_for_file, start benchmarking\n");
    // Benchmarking
    int r;
    for (r = 0; r < repeats; r++) {
	//printf("DEBUG: benchmark round %d\n", r);
	bool success = false;

        // Start time
        double start_time = omp_get_wtime();
        
        // Invert matrix
        double** mat_inv = invert_matrix(n, n, mat_cp, &success);
        
        // End time
        double end_time = omp_get_wtime();
        
        double elapsed = end_time - start_time;
	printf("rep %d for mat %d: start %f, end %f\n", r, n, start_time, end_time);
        //save_results("./HPC.ParallelMatrixInversion/results.csv", n, num_threads, r + 1, elapsed);        
        save_results(save_location, n, num_threads, r + 1, elapsed);

        // Check the inverse for fun
        if (!check_inverse(n, n, mat, mat_inv)) {
            printf("Failed to find inverse for %s\n", fname);
        }
	free_mat(mat_inv, n);
    }

    free_mat(mat, n);
    free_mat(mat_cp,  n);
    return true;
}

void save_results(const char* fname, int size, int threads, int run, double time) {
    FILE* file = fopen(fname, "a");
    if (!file) {
        fprintf(stderr, "Error opening the result file\n");
        return;
    }
    fprintf(file, "%d,%d,%d,%f\n", size, threads, run, time);
    fclose(file);
    //printf("DEBUG: %d, %d, %d, %f saved to file %s\n", size, threads, run, time, fname);
}


int run_test_matrices() {
    print_working_dir();

    const char *directory_path = "./HPC.ParallelMatrixInversion/test_matrices/";
    //const char *directory_path = "test_matrices/"; // Debugging path
    DIR *dir = opendir(directory_path);
    if (!dir) {
        perror("opendir");
        return 1;
    }

    int failed_tests = 0;

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        /* Skip directories (there shouldn't be any) */
        if (entry->d_type == DT_REG) {
            if (!invert_mat_from_file(directory_path, entry->d_name)) {
                failed_tests += 1;
            }
        }
    }
    closedir(dir);

    return failed_tests;
}

/* A method to test invert_matrix with specific set of test matrices 
 * The test matrices are either invertible or singular, which is determined by their
 * filename: the format is mat_nxn_c, where n is the dimension and c is either i (invertible)
 *  or s (singular)
*/
bool invert_mat_from_file(const char* dir, const char* fname) {

    /* Create the target file path  */
    char* full_path = malloc(strlen(dir) + strlen(fname) + 1);

    if (!full_path) {
        perror("malloc (concatenate file path)");
        return false;
    }

    strcpy(full_path, dir);
    strcat(full_path, fname);

    /* Extract dimensions and type of matrix (invertible or singular) from the filename */
    printf("Reading matrix from file %s ", fname);
    int nrow, ncol;
    char kind;

    if (sscanf(fname, "mat_%dx%d_%c", &nrow, &ncol, &kind) != 3) {
        printf("Invalid filename format: %s\n", fname);
	free(full_path);
        return false;
    }

    // printf("Reading %dx%d matrix from %s ", nrow, ncol, full_path);
    printf("(%s)\n", (kind == 'i') ? "invertible" : "singular");

    /* Open the file */
    FILE *fp = fopen(full_path, "r");
    if (!fp) {
        perror("Error opening file");
        print_working_dir();
	free(full_path);
        return false;
    }

    /* Scan the file and read the matrix */
    int i, j;
    double** mat = allocate_mat(nrow, ncol);
    if (!mat) {
	printf("Error in invert_mat_from_file: memory allocatoin (mat) failed\n");
	free(full_path);
	return false;
    }

    for (i = 0; i < nrow; ++i) {
        for (j = 0; j < ncol; ++j) {
            if (fscanf(fp, "%lf", &mat[i][j]) != 1) {
                printf("Error reading matrix value at [%d][%d] in file: %s\n", i, j, fname);
                fclose(fp);
		free(full_path);
		free_mat(mat, nrow);
                return false;
            }
        }
    }

    /* Close the file after successfully reading the matrix */
    fclose(fp);
    free(full_path);

    //printf("DEBUG: Matrix read from file successfully");

    /* Create a copy of the input matrix before inverting it */
    //double mat_cp[nrow][ncol];
    double** mat_cp = allocate_mat(nrow, ncol);
    if (!mat_cp) {
	printf("Error in invert_mat_from_file: memory allocation (mat_cp) failed\n");
	return false;
    }

    copy_mat(nrow, ncol, mat, mat_cp);

    /* Invert the copy of original matrix */
    //double mat_inv[nrow][ncol];
    //double** mat_inv = allocate_mat(nrow, ncol);
    bool res = false;
    double** mat_inv = invert_matrix(nrow, ncol, mat_cp, &res);

    if (res && kind == 'i') {
        if (check_inverse(nrow, ncol, mat, mat_inv)) {
            printf("Found correct inverse for %s\n", fname);
        } else {
            printf("ERROR: Didn't find inverse for invertible matrix at %s\n", fname);
        }
	free_mat(mat_inv, nrow); // Free the inverse only if it was found
    } else if (!res && kind == 'i') {
        printf("ERROR: Failed to invert invertible matrix at %s\n", fname);
    }
    
    //printf("DEBUG: Freeing the matrices ");
    free_mat(mat, nrow);
    free_mat(mat_cp, nrow);
    //free_mat(mat_inv, nrow);
    printf("\n\n");
    return true;
}

/* Function to test that the found inverse is the actual inverse. Used as a test for generated
 * matrices */
bool check_inverse(int nrow, int ncol, double** m1, double** m2) {
    bool success = true;
    //double mat_res[nrow][ncol];
    double** mat_res = allocate_mat(nrow, ncol);
    if (!mat_res) {
	printf("Error in check_inverse: memory allocation (mat_res) failed\n");
	return false;
    }

    int i, j, k;

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            mat_res[i][j] = 0;

            for (k = 0; k < nrow; k++) {
                mat_res[i][j] += m1[i][k] * m2[k][j];
            }

            /* Check if the result is identity (will fail if there are nans or +-inf) */
            if (i == j && fabs(1 - mat_res[i][j]) > 1e-9) {
                success = false;
            }

            if (i != j && fabs(mat_res[i][j]) > 1e-9) {
                success = false;
            }
        }
    }

//    if (!success) {
//        printf("A*A^-1:\n");
//        print_mat(nrow, nrow, mat_res);
//    }
    free_mat(mat_res, nrow);
    return success;
}

