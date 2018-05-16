/* ---------------------main.c--------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "gaussian.h"


int main()
{
	double *originA, *A, *pA;  // matrix   A for single thread; pA for OpenMP
	double *x, *y, *xnew, *ynew;   // vectors x,y for single thread; xnew, ynew for OpenMP
	double *b, *pb, *originb;  // vectors  b for single thread; pb for OpenMP
	double start, end; //time conuter
	int r, c; // row and column 
	int n; //number of equations
	int thread_count, correct_count, correct_count_omp;
	int result;

	printf("Please input the number of equations:");
	scanf("%d", &n);

	r = n;
	c = n;

	originA = (double*)malloc(r * c * sizeof(double)); //memory for original Matrix
	if (originA == NULL)
		return -1;

	A = (double*)malloc(r*c*sizeof(double));	//memory for Matrix
	if (A == NULL)
		return -1;

	pA = (double*)malloc(r*c*sizeof(double));	//memory for Matrix
	if (pA == NULL)
		return -1;


	originb = (double*)malloc(r * sizeof(double)); //memory for original vector
	if (originb == NULL)
		return -1;

	b = (double*)malloc(r * sizeof(double)); //memory for vector
	if (b == NULL)
		return -1;
	
	pb = (double*)malloc(r * sizeof(double)); //memory for vector
	if (pb == NULL)
		return -1;	



	x = (double *)malloc(r * sizeof(double));
	if (x == NULL)
		return -1;

	y = (double*)malloc(r * sizeof(double)); //memory for substitution vector
	if (y == NULL)
		return -1;

	xnew = (double *)malloc(r * sizeof(double));
	if (xnew == NULL)
		return -1;

	ynew = (double*)malloc(r * sizeof(double)); //memory for substitution vector
	if (ynew == NULL)
		return -1;


	/* Initialize the matrix X and vector v*/
	srand((unsigned)time(NULL));

	//inintialize matrix and vector
	Mat_Init(r, c, originA);
	Vec_Init(r, originb);

	//copy matrix from A to originA and pA;
	for (int i = 0; i < n * n; i++)
	{
		A[i] = originA[i];
		pA[i] = originA[i];
	}
	//copy vector from B to originA and pb;
	for (int i = 0; i < n; i++)
	{
		b[i] = originb[i];
		pb[i] = originb[i];
	}

	printf("\noriginal Matrix:");
//	Mat_Show(n, n, originA);
//  Vec_Show(n, originb);


	printf("\n-----------------------------------------------\nCalculation implements the single thread: ");
	//printf("\n--After elimination------------------\n");
	/* calculate and record its running time */
	start = omp_get_wtime();
	gauss_elimination(A, n, b, x, y);
	end = omp_get_wtime();
	
	printf("Single thread takes %f millisecond.\n", (end - start) * 1000);

	//do verification
	correct_count = verification(x, originA, originb, n);
	printf("%d ", correct_count);
	printf("The single thread calculation result is: ");
	if (correct_count == n)
	{
		printf("right\n!");
	}
	else
	{
		printf("wrong\n!");
	}
	


	//using openMP 	
	printf("\n-----------------------------------------------\nCalculation implements the OpenMP: ");
	
		//try do many times for one set of equations
		printf("\nPlease input the number of threads: ");
		scanf("%d", &thread_count);
		/* calculate and record its running time */
		start = omp_get_wtime();
		gauss_elimination_omp(pA, n, pb, xnew, ynew, thread_count);
		end = omp_get_wtime();
		
		
		printf("The calculation cost %f milliseconds with openmp methord\n", (end - start) * 1000);
				//do verification
		correct_count_omp = verification(xnew, originA, originb, n);
		printf("%d ", correct_count_omp);
		printf("The multi-thread calculation result is: ");
		if (correct_count_omp == n)
		{
			printf("The result is right\n!");
		}
		else
		{
			printf("The result is wrong\n!");
		}
	

	free(A);
	free(b);
	free(x);
	free(y);
	free(pA);
	free(pb);
	free(xnew);
	free(ynew);
	return(0);
}

//function fo verifycation of the result
int verification(double *x, double *A, double *originb, int n)
{
	double *verifyB;
	int correctRNo = 0;

	verifyB = (double *)malloc(n * sizeof(double));
	if (verifyB == NULL)
		return -1;

	//caculate the verifid b with the calculated result of x[i]
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			verifyB[i] += A[i * n + j] * x[j];
		}
	}

	//increase the correct NO. count  
	for (int i = 0; i < n; i++)
	{
		if (fabs(verifyB[i] - originb[i]) <= 0.000001)
		{
			correctRNo = correctRNo + 1;
			printf("i = %d: %f %f\n", i, verifyB[i], originb[i]);
		}
		else{
			printf("i = %d: %f %f\n", i, verifyB[i], originb[i]);
		}
	}

	free(verifyB);
	return correctRNo;
}
