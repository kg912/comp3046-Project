/* ---------------------gaussian.h--------------------- */
#ifndef GAUSSIAN_H
#define GAUSSIAN_H

/* initialize a matrix x */
void Mat_Init(int row, int col, double *X);

/* display a matrix */
void Mat_Show(int row, int col, double *A);

/* initialize a vector x */
void Vec_Init(int size, double *x);

/* display a vector */
void Vec_Show(int size, double *b);

/* gaussian elimination function using single thread */
void gauss_elimination(double *A, int n, double *b, double *x, double *y);

/* gaussian elimination function using parallel threads r */
void gauss_elimination_omp(double *pA, int n, double *pb, double *x, double *y, int thread_count);

#endif










