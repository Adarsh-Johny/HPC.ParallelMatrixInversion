#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

// Configuration for number of matrices for each size
// Format: { {matrix_size, count}, ... }
int matrix_config[][2] = {
    {2, 1},
    {3, 1},
    {4, 1},
    {5, 1},
    {6, 1},
    {7, 1},
    {8, 1},
    {9, 1},
    {10, 1}};
int config_size = sizeof(matrix_config) / sizeof(matrix_config[0]);

// Function to calculate the determinant of a matrix (recursive for general case)
double determinant(double **matrix, int n)
{
    if (n == 1)
    {
        return matrix[0][0];
    }

    if (n == 2)
    {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }

    double det = 0;
    for (int p = 0; p < n; p++)
    {
        double **submatrix = (double **)malloc((n - 1) * sizeof(double *));
        for (int i = 0; i < n - 1; i++)
        {
            submatrix[i] = (double *)malloc((n - 1) * sizeof(double));
        }

        for (int i = 1; i < n; i++)
        {
            int col_index = 0;
            for (int j = 0; j < n; j++)
            {
                if (j == p)
                    continue;
                submatrix[i - 1][col_index++] = matrix[i][j];
            }
        }

        det += (p % 2 == 0 ? 1 : -1) * matrix[0][p] * determinant(submatrix, n - 1);

        for (int i = 0; i < n - 1; i++)
        {
            free(submatrix[i]);
        }
        free(submatrix);
    }

    return det;
}

// Function to save the matrix to a file
void save_matrix_to_file(double **matrix, int size, int index)
{
    char folder_name[] = "../test_matrices";
    mkdir(folder_name, 0777); // Create the folder if it doesn't exist

    char file_name[100];
    sprintf(file_name, "%s/matrix_%dx%d_%02d.txt", folder_name, size, size, index);

    FILE *file = fopen(file_name, "w");
    if (file == NULL)
    {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            fprintf(file, "%.2f ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

// Function to generate a random matrix and ensure it is invertible
void generate_invertible_matrix(int size, int min_value, int max_value, int num_matrices)
{
    double **matrix;
    matrix = (double **)malloc(size * sizeof(double *));
    for (int i = 0; i < size; i++)
    {
        matrix[i] = (double *)malloc(size * sizeof(double));
    }

    for (int idx = 1; idx <= num_matrices; idx++)
    {
        while (1)
        {
            // Fill the matrix with random values
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    matrix[i][j] = (rand() % (max_value - min_value + 1)) + min_value;
                }
            }

            // Check if the matrix is invertible
            if (fabs(determinant(matrix, size)) > 1e-6)
            {
                break;
            }
        }

        // Save the matrix to a file
        save_matrix_to_file(matrix, size, idx);
    }

    // Free allocated memory
    for (int i = 0; i < size; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

int main()
{
    srand(time(0)); // Seed for random number generation

    for (int i = 0; i < config_size; i++)
    {
        int size = matrix_config[i][0];
        int count = matrix_config[i][1];
        generate_invertible_matrix(size, 0, 100, count);
    }
    return 0;
}
