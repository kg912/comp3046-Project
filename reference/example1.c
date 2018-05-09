#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Windows.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include <omp.h> 

int checkans(int a, double *oA, double *x, double *b)
{
	int i, r;
	int result = 0;
	double *re;

	re = (double *)malloc(a * sizeof(double));
	if (re == NULL)
		return -1;

	for ( r=0 ; r < a; r++)
	{
		re[r] = 0;
	}

	for (int p = 0; p < a; p++)
	{
		for (int q = 0; q < a; q++)
		{
			re[p] += oA[p * a + q] * x[q];
		}
	}

	for (i = 0; i < a; i++)
	{
		if (fabs(b[i] - re[i]) <= 0.0001)
		{
			result = result + 1;
		}
	}

	free(re);
	return (result);
}

int checkansO(int a, double *oA, double *x, double *b, int thread_count)
{
	int i, r, n,q,p;
	int result = 0, myresult = 0;
	double *re;

	re = (double *)malloc(a * sizeof(double));
	if (re == NULL)
		return -1; 
	for (r = 0; r < a; r++)
	{
		re[r] = 0;
	}

//first	
#  pragma omp parallel for num_threads(thread_count) \
	private(p,q)
	for (p = 0; p < a; p++)
	{
		for (q = 0; q < a; q++)
		{
			re[p] += oA[p * a + q] * x[q];
		}
	}
//second
#  pragma omp parallel for num_threads(thread_count) \
        reduction(+:result) private(i)
	for (i = 0; i < a; i++)
	{
		if (fabs(b[i] - re[i]) <= 0.0001)
		{
			result = result + 1;
		}
	}

	free(re);
	return (result);
}

int main()
{
	double *A, *B, *y, *x, *oA, *oB; /* A and B */
	DWORD start, end;
	int a, b;                                /*size of A and B*/
	int i, j, k, n;
	int result,index;

	int resultO;
	int thread_count;
	double startO, endO;
	double *AO, *BO, *yO, *xO, *oAO, *oBO;
	double *tmp, *tmpO; /* A and B */
	double max, maxO;



	//generate A and B with input size
	srand((unsigned)time(NULL));
	printf("Please input a number: \n");
	scanf("%d", &a);

	A = (double *)malloc(a * a * sizeof(double));
	if (A == NULL)
		return -1;

	oA = (double *)malloc(a * a * sizeof(double));
	if (oA == NULL)
		return -1;

	B = (double *)malloc(a * sizeof(double));
	if (B == NULL)
		return -1;

	oB = (double *)malloc(a * sizeof(double));
	if (oB == NULL)
		return -1;

	x = (double *)malloc(a * sizeof(double));
	if (x == NULL)
		return -1;

	y = (double *)malloc(a * sizeof(double));
	if (y == NULL)
		return -1;

	oAO = (double *)malloc(a * a * sizeof(double));
	if (oAO == NULL)
		return -1;

	oBO = (double *)malloc(a * sizeof(double));
	if (oBO == NULL)
		return -1;

	xO = (double *)malloc(a * sizeof(double));
	if (xO == NULL)
		return -1;

	yO = (double *)malloc(a * sizeof(double));
	if (yO == NULL)
		return -1;

	AO = (double *)malloc(a * a * sizeof(double));
	if (AO == NULL)
		return -1;

	BO = (double *)malloc(a * sizeof(double));
	if (BO == NULL)
		return -1;

	tmp = (double *)malloc(a * sizeof(double));
	if (tmp == NULL)
		return -1;

	tmpO = (double *)malloc(a * sizeof(double));
	if (tmpO == NULL)
		return -1;

	for (int w = 0; w < a; w++)
	{
		x[w] = 0;
		y[w] = 0;
		xO[w] = 0;
		yO[w] = 0;
	}

	for (int e = 0; e < a * a; e++)
	{
		oA[e] = 0;
		oAO[e] = 0;
	}

	Mat_Init(a, a, A);
	Vec_Init(a, B);

	for (int q = 0; q < a * a; q++)
	{
		oA[q] = A[q];
		oAO[q] = A[q];
		AO[q] = A[q];
	}

	for (int q = 0; q < a; q++)
	{
		oB[q] = B[q];
		BO[q] = B[q];
		oBO[q] = B[q];
	}

	// without using openmp
	start = GetTickCount();

	//partial pivoting
	index = 0;
	max = A[0];
	for (k = 0; k < a; k++)
	{
		// printf("The max is: %f\n", max);
		// printf("The A is: %f\n", A[k * a]);
		if (A[k * a] > max)
		{
			index = k;
			max = A[k * a];
			// printf("The index changed to: %d\n", index);
			// printf("The max changed to: %f\n", A[k * a]);
		}
	}

	//printf("max index: %d\n", index);
	//printf("max: %f\n", A[index * a]);

	for (k = 0; k < a; k++)
	{
		tmp[k] = A[k];
		A[k] = A[index * a + k];
		A[index * a + k] = tmp[k];
	}

	double tmp2 = 0.0;
	tmp2 = B[0];
	B[0] = B[index];
	B[index] = tmp2;

	for (k = 0; k < a; k++)
	{
		for (j = k + 1; j < a; j++)
		{
			A[k * a + j] = A[k * a + j] / A[k + k * a];
		}
		y[k] = B[k] / A[k + k * a];

		A[k + k * a] = 1;
		for (i = k + 1; i < a; i++)
		{
			for (j = k + 1; j < a; j++)
			{
				A[i * a + j] = A[i * a + j] - A[i * a + k] * A[k * a + j];
			}
			B[i] = B[i] - A[i * a + k] * y[k];
			A[i * a + k] = 0;
		}
	}
	for (k = a - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (i = k - 1; i >= 0; i--)
		{
			y[i] = y[i] - x[k] * A[i * a + k];
		}
	}
    result = checkans(a, oA, x, oB);
    printf("%d ", result);
    if (result == a)
    {
        printf("The CPU ans is right\n!");
    }
    else
    {
        printf("The CPU ans is wrong\n!");
    }
	end = GetTickCount();
	printf("It takes %d milliseconds for no openmp\n", (end - start));


	//using openmp 
	printf("Please input the number of threads: ");
	scanf("%d", &thread_count);

	startO = omp_get_wtime();

	//partial pivoting
	index = 0;
	max = A[0];
	for (k = 0; k < a; k++)
	{
		// printf("The max is: %f\n", max);
		// printf("The A is: %f\n", A[k * a]);
		if (A[k * a] > max)
		{
			index = k;
			max = A[k * a];
			// printf("The index changed to: %d\n", index);
			// printf("The max changed to: %f\n", A[k * a]);
		}
	}

	//printf("max index: %d\n", index);
	//printf("max: %f\n", A[index * a]);

	for (k = 0; k < a; k++)
	{
		tmp[k] = A[k];
		A[k] = A[index * a + k];
		A[index * a + k] = tmp[k];
	}

	tmp2 = 0.0;
	tmp2 = B[0];
	B[0] = B[index];
	B[index] = tmp2;

	for (k = 0; k < a; k++)
	{
//the second
#  pragma omp parallel for num_threads(thread_count)
		for (i = 0; i < thread_count; i++){
			for (j =(k + 1)/thread_count; j < a/thread_count; j++)
			{
				AO[k * a + j] = AO[k * a + j] / AO[k + k * a];
			}
		}
		yO[k] = BO[k] / AO[k + k * a];
		AO[k + k * a] = 1;
//the third
#  pragma omp parallel for num_threads(thread_count) \
        default(none) private(i,j) shared(a,AO,BO,yO,k)
		for (i = k + 1; i < a; i++)
		{
			for (j = k + 1; j < a; j++)
			{
				AO[i * a + j] = AO[i * a + j] - AO[i * a + k] * AO[k * a + j];
			}
			BO[i] = BO[i] - AO[i * a + k] * yO[k];
			AO[i * a + k] = 0;
		}
	}

	for (k = a - 1; k >= 0; k--)
	{
		xO[k] = yO[k];
//forth
#  pragma omp parallel for num_threads(thread_count) \
	default(none) private(i) shared(yO,xO,AO,a,k)
		for (i = k - 1; i >= 0; i--)
		{
			yO[i] = yO[i] - xO[k] * AO[i*a + k];
		}
	}

	resultO = checkansO(a, oAO, xO, oBO, thread_count);
	//resultO = checkans(a, oAO, xO, oBO);
    printf("%d ", resultO);
    if (resultO == a)
    {
        printf("The CPU using openmp ans is right\n!");
    }
    else
    {
        printf("The CPU using openmp ans is wrong\n!");
    }
	endO = omp_get_wtime();
	printf("It takes %f milliseconds for using openmp\n", (endO - startO) * 1000);
	
	free(A);
	free(oA);
	free(B);
	free(x);
	free(oB);
	free(y);
	free(yO);
	free(AO);
	free(BO);
	free(oAO);
	free(oBO);
	free(xO);

	return (0);
}

