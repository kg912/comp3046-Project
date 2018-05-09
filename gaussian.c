// =============================================================================
#include <stdio.h>
#include <stdlib.h>
#include "gaussian.h"
// =============================================================================


/* initialize a matrix x. You need to call srand() in advance. */
void Mat_Init(int row, int col, double *X)
{
	int i, size;

	size = row * col;
	for (i = 0; i < size; i++)
		X[i] = (double)((float)rand()) / ((float)RAND_MAX) * 100;
}



/* display a matrix */
void Mat_Show(int row, int col, double *X)
{
	int i, j;
	printf("row = %d col = %d\n", row, col);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++)
			printf("%.4f ", X[col * i + j]);
		printf("\n");
	}
}


/* initialize a vector x */
void Vec_Init(int size, double *x)
{
	int i;

	for (i = 0; i < size; i++)
		x[i] = (double)((float)rand()) / ((float)RAND_MAX) * 100;
}

/* display a vector */
void Vec_Show(int size, double *x)
{
	int i;

	printf("\nvector size = %d\n", size);
	for (i = 0; i < size; i++)
		printf("%.4f ", x[i]);
	printf("\n");
}


//free a matrix
void free_matrix(double **m, int r, int c)
{
	int i;

	if (m) {
		for (i = 0; i < r; i++) {
			free(m[i]);
		}
		free(m);
	}
}



void gauss_elimination(double *A, int n, double *b, double *x, double *y)
{
	int i, j, k;
	double temp;
	//gauss elimination parts
	for (k = 0; k < n - 1; k++) {
		for (j = k + 1; j < n - 1; j++) {
			temp = A[k * n + j] / A[k * n + k];
		}
		y[k] = b[k] / A[k *n + k];
		A[k *n + k] = 1;

		for (i = k + 1; i <= n - 1; i++) {
			for (j = k + 1; j < n - 1; j++){
				A[i*n + j] = A[i*n + j] - A[i* n + k] * temp;
			};
			b[i] = b[i] - A[i * n + k] * y[k];
			A[i * n + k] = 0;
		}
	}
	printf("\nMatrix after elimination: ");
	
	Mat_Show(n, n, A);
	printf("\nvectors after elimination:");
	Vec_Show(n, b);
	
	//Back-substitution parts
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (i = k - 1; i >= 0; i--)
		{
			y[i] = y[i] - x[k] * A[i * n + k];
		}
	}
	printf("\n-----------------------------------------------\nResults: ");
	Vec_Show(n, x);
}
			
	
