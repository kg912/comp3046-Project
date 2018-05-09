/* matrix.h */

#ifndef GAUSSIAN_H
#define GAUSSIAN_H


// =============================================================================
/* initialize a matrix x */
void Mat_Init(int row, int col, double *X);

/* display a matrix */
void Mat_Show(int row, int col, double *A);

/* initialize a vector x */
void Vec_Init(int size, double *x);

/* display a vector */
void Vec_Show(int size, double *b);

// =============================================================================
double *calculate_solution(double **m, int r);
// =============================================================================
void fill_matrix(double **m, int r, int c, char *file);
// =============================================================================
void print_matrix(double **m, int r, int c, char *file);
// =============================================================================
void free_matrix(double **m, int r, int c);
// =============================================================================
void gauss_elimination(double *A, int n, double *b, double *x, double *y);
// =============================================================================


#endif
