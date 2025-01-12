#ifndef FILE_READER_H
#define FILE_READER_H

#include <stdbool.h> /* bool */

/* Reads a matrix from a file.
 *
 * dir: Directory path to the file.
 * fname: Name of the file.
 * nrow: Pointer to store the number of rows of the matrix.
 * ncol: Pointer to store the number of columns of the matrix.
 * mat: Pointer to store the dynamically allocated matrix.
 *
 * Returns true on success, false on failure.
 */
bool read_matrix_from_file(const char *filepath, int *nrow, int *ncol, double ***mat);

#endif /* FILE_READER_H */
