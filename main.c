#include <stdio.h>
#include <windows.h>
#include <time.h>
#include <omp.h>
#include "gaussian.h"


int main()
{
	double *A, *originA, *tmp; /* matrix */
	double *y, *b, *originb; /* vectors */
	double start, end;
	int r, c;
	int n; //number of equations
	int thread_count;

	printf("Please input the number of equations:");
	scanf("%d", &n);

	r = n;
	c = n;

	A = (double*)malloc(r*c*sizeof(double));
	if (A == NULL)
		return -1;

	originA = (double*)malloc(r*c*sizeof(double));
	if (originA == NULL)
		return -1;

	b = (double*)malloc(r*sizeof(double));
	if (b == NULL)
		return -1;

	originb = (double*)malloc(r*sizeof(double));
	if (originb == NULL)
		return -1;

	y = (double*)malloc(r*sizeof(double));
	if (y == NULL)
		return -1;

	originA = A;

	/* Initialize the matrix X and vector v*/
	srand((unsigned)time(NULL));

	originA = A;
	originb = b;

	Mat_Init(r, c, A);
	Vec_Init(r, b);

	printf("\n Original Matrix: ");
	Mat_Show( r,  c,  A);
	printf("\n Original vectors: ");
	Vec_Show(n , b);
	printf("\n After elimination: ");
	gauss_elimination(A, n, b, y);
	Mat_Show(n, n, A);
	Vec_Show(n, b);
	/* calculate Xv and record its running time */
	start = omp_get_wtime();

	end = omp_get_wtime();
}