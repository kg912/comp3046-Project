#include <stdio.h>
#include <windows.h>
#include <time.h>
#include <omp.h>
#include "gaussian.h"


int main()
{
	double *A, *originA, *tmp; /* matrix */
	double *x, *y, *b, *originb; /* vectors */
	double start, end; //time conuter
	int r, c; // row and column 
	int n; //number of equations
	int thread_count; 

	printf("Please input the number of equations:");
	scanf("%d", &n);

	r = n;
	c = n;

	A = (double*)malloc(r*c*sizeof(double));	//memory for Matrix
	if (A == NULL)
		return -1;

	originA = (double*)malloc(r*c*sizeof(double)); //memory for original Matrix
	if (originA == NULL)
		return -1;

	b = (double*)malloc(r*sizeof(double)); //memory for vector
	if (b == NULL)
		return -1;

	originb = (double*)malloc(r*sizeof(double)); //memory for original vector
	if (originb == NULL)
		return -1;

	x = (double *)malloc(r * sizeof(double));
	if (x == NULL)
		return -1;

	y = (double*)malloc(r*sizeof(double)); //memory for substitution vector
	if (y == NULL)
		return -1;

	originA = A;

	/* Initialize the matrix X and vector v*/
	srand((unsigned)time(NULL));

	originA = A;
	originb = b;

	//inintialize matrix and vector
	Mat_Init(r, c, A);
	Vec_Init(r, b);

	printf("\n-----------------------------------------------\nOriginal Matrix: ");
	Mat_Show( r,  c,  A);
	printf("\n-----------------------------------------------\nOriginal vectors:");
	Vec_Show(n , b);
	printf("\n----------------After elimination------------------\n");
	gauss_elimination(A, n, b, x, y);
	
	/* calculate Xv and record its running time */
	start = omp_get_wtime();

	end = omp_get_wtime();
}
