#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"otraz.h"


int main (int argc, char** argv)
{
	if (!(argc == 4 || argc == 5))
	{
		printf ("USAGE: ./a.out	 size_matr  number_of_number  EPS  file_name\n or\n./a.out  size_matr  number_of_number  EPS\n");
		return -1;
	}
	int n, k;
	double *a;
	double eps;
	n = atoi(argv[1]);
	
	k = atoi(argv[2]);
	eps = atof(argv[3]);
	
	if (n <= 0 || k <= 0 || eps < 0){
		printf("Incorect data\n");
		return -11;
	}
		
	a = new double [n* n];
	
	if (a == NULL)
	{
		printf("ERROR of memory\n");
		return -2;
	}
	
	if (argc == 4)
	{
		read_func(a,n);
	}
	else
	{
		if (read_file(a,n,argv) < 0)
		{
			printf("ERROR in the reading function\n");
			delete[] a;
			return -3;
		}
			
	}
		
	//printf("\n");
	
	double t1 = clock ();
	
	mirror (a,n);
	printf ("TIME_mirror: %f\n", (clock() - t1)/CLOCKS_PER_SEC);
    //print_mat (a,n);

	double *A = new double [n];
	if (A == NULL)
	{
		printf ("ERROR of memory\n");
		return -3;	
	}
	
	double *B = new double [n];
	if (B == NULL)
	{
		printf ("ERROR of memory\n");
		return -4;	
	}
	
	for (int i = 0; i < n; i++)
	{
		A[i] = a[i * n + i];
	}
	
	for (int i = 0; i < (n - 1); i++)
	{
		B[i] = a[i * n + (i + 1)];
	}	
	delete[] a;
	
	/*
	double help;	
	scanf ("%lf", &help);
	n_l (A, B, n, help);	
	*/
	t1 = clock();		
    if (solution (A, B, n, k, eps) < 0)
	{
		printf ("ERROR\n");
		delete[] A;
		delete[] B;
		
		return -5;
    }
    /*double l;
    scanf ("%lf", &l);
    n_l (A, B, n,l);
    */
	printf ("TIME_bissection: %f\n", (clock() - t1)/CLOCKS_PER_SEC);
		
	/*
	for (int i = 0; i < n; i++)
	{
		printf ("%f ", A[i]);
	}
	printf ("\n");
	for (int i = 0; i < n - 1; i++)
	{
		printf ("%f ", B[i]);
	}
	printf ("\n");
	*/
	
		
	delete[] A;
	delete[] B;				
	return 0;
}
