#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "gaussian.h"

int main()
{
	double  *originA, *A, *B; /* matrix */
	double *x, *y, *b, *originb, *r1, *r2; /* vectors */

	double start, end, start1, end1; //time conuter
	int r, c; // row and column 
	int n; //number of equations
	int thread_count;
	int result;

	printf("Please input the number of equations:");
	scanf("%d", &n);

	r = n;
	c = n;

	A = (double*)malloc(r*c*sizeof(double));	//memory for Matrix
	if (A == NULL)
		return -1;

	originA = (double*)malloc(r * c * sizeof(double)); //memory for original Matrix
	if (originA == NULL)
		return -1;

	B = (double*)malloc(r*c*sizeof(double));	//memory for Matrix
	if (A == NULL)
		return -1;

	b = (double*)malloc(r * sizeof(double)); //memory for vector
	if (b == NULL)
		return -1;

	originb = (double*)malloc(r * sizeof(double)); //memory for original vector
	if (originb == NULL)
		return -1;



	x = (double *)malloc(r * sizeof(double));
	if (x == NULL)
		return -1;

	y = (double*)malloc(r * sizeof(double)); //memory for substitution vector
	if (y == NULL)
		return -1;

	r1 = (double *)malloc(r * sizeof(double)); //memory for result of single thread calculation
	if (r1 == NULL)
		return -1;

	r2 = (double *)malloc(r * sizeof(double)); //memory for result OpenMp calculation
	if (r1 == NULL)
		return -1;

	/* Initialize the matrix X and vector v*/
	srand((unsigned)time(NULL));

	//inintialize matrix and vector
	Mat_Init(r, c, A);
	Vec_Init(r, b);

	//copy matrix from A to originA;
	for (int i = 0; i < n * n; i++)
	{
		originA[i] = A[i];
		B[i] = A[i];
	}
	//copy vector from B to originA;
	for (int i = 0; i < n; i++)
	{
		originb[i] = b[i];
	}

	printf("\n-----------------------------------------------\noriginal Matrix:");
	Mat_Show(n, n, A);
	Vec_Show(n, b);


	printf("\n-----------------------------------------------\nCalculation implements the single thread: ");
	printf("\n--After elimination------------------\n");
	/* calculate and record its running time */
	start = omp_get_wtime();
	gauss_elimination(A, n, b, x, y);
	end = omp_get_wtime();
	printf("Single thread takes %f millisecond.\n", (end - start) * 1000);
	for (int i = 0; i < n; i++)
	{
		r1[i] = x[i];
	}
	//do verification
	int correctRNo = verification(r1, r2);
	printf("%d ", correctRNo);
	if (correctRNo == n)
	{
		printf("The result is right\n!");
	}
	else
	{
		printf("The result is wrong\n!");
	}
	
	//using openMP 	
	printf("\n-----------------------------------------------\nCalculation implements the OpenMP: ");
	printf("\nIt will loop 10 times to test withdiffernet thread counts");

	//try do many times for one set of equations
	for (int m = 0; m < 10; m = m + 1){
	
		printf("\nPlease input the number of threads: ");
		scanf("%d", &thread_count);
		/* calculate and record its running time */

		start = omp_get_wtime();
		gauss_elimination_omp(A, n, b, x, y, thread_count);
		end = omp_get_wtime();
		printf("The calculation cost %f milliseconds with openmp methord\n", (end - start) * 1000);
		for (int i = 0; i < n ; i++)
		{
			r2[i] = x[i];
		}

		//do verification
		int correctRNo = verification(r1, r2 );
		printf("%d ", correctRNo);
		if (correctRNo == n)
		{
			printf("The result is right\n!");
		}
		else
		{
			printf("The result is wrong\n!");
		}	
	}

	free(A);	
	free(b);
	free(x);
	free(y);
	return(0);
}

//function fo verifycation of the result
int verification(double *x, double *y, int n)
{
	int correctRNo = 0;

	//caculate the verifid b with the calculated result of x[i]
	for (int i = 0; i < n; i++)
	{		
		if (fabs(x[i] - y[i]) <= 0.0000001)//verifycation
		{
			correctRNo = correctRNo + 1;   //increase the correct NO. count  
		}
	}

	return correctRNo;
}
