#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"otraz.h"


int solution (double *main, double * sub_main, int n, int k, double eps)
{
	double a_0 = 0, b_0 = 0;
	
	double *A = new double [n];
	if (A == NULL){
		printf ("Error of memory\n");
		return -1;
	}
	
	double *B = new double [n];
	if (B == NULL){
		printf("Error of memory\n");
		delete[] A;
		return -2;
	}
	
    int i = 0;
    while (i < 50){
		a_0 -= 100;
		b_0 += 100;
		copy_array (main, A, n);
		copy_array (sub_main, B, n - 1);
		
		if (n_l (A, B, n, a_0) < k)
		{
			copy_array (main, A, n);
			copy_array (sub_main, B, n - 1);
			if (n_l (A, B, n, b_0) >= k)
				break;
		}
		i++;
	
	}
    printf("HERE\n\n");

    if (i >= 50){
		printf("Can't find sector\n");
		delete[] A;
		delete[] B;
		return -3;
	}
	
	double c;
		
	while (b_0 - a_0 > eps)
	{
		c = a_0 + b_0;
		c /= 2;
		copy_array (main, A, n);
		copy_array (sub_main, B, n - 1);
			
		if (n_l(A, B, n, c) < k)
			a_0 = c;
		else
			b_0 = c;
//		printf ("%f\n", c);
	}
	printf ("answer: %f\n", (a_0 + b_0) / 2.0);

	delete[] A;
	delete[] B;
    return 1;
}

void copy_array (const double *from, double *a,int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = from[i];
	}
}


int n_l (double *main, double *sub_main, int n, double l)
{
	for(int i = 0; i < n; i++){
		main[i] -= l;
	}
	
	change_matr (main, sub_main, n);	
    long double x = main[0];
    long double y = 1;
    long double a, b;
    long double u, v;
	int m;
	if (x * y < 0)
	{
		m = 1;
	}
	else
		m = 0;
    long double gamma;
	for (int k = 1; k < n; k++)
	{
		a = main[k];
        b = sub_main[k - 1];
		if (b < EPS)
			b = 0.0;
        gamma = (1.0 / EPS) * fmax (fabs (x), fabs (b* (b * y)));
		u = gamma * (a * x - b * b * y);
		v = gamma * x;
        if (!((u > 0 && x > 0 ) || (u < 0 && x < 0)))
			m++;
		x = u;
		y = v; 
	}
    //printf ("%d\n", m);
	return m;
}

void change_matr (double *A, double *B, int n){
	double max1 = fabs (A[0]), max2 = fabs (B[0]);
	for (int i = 0; i < n - 1; i++)
	{
		if(max1 < fabs (A[i]))
			max1 = fabs (A[i]);
		if(max2 < fabs (B[i]))
			max2 = fabs (B[i]);
	}
	
	if (max1 < fabs (A[n-1]))
		max1 = fabs (A[n-1]);
		
	if (max2 > max1)
		max1 = max2;

	max1 *= 4;
	//printf ("hello %f\n", max1);
	for (int i = 0; i < (n - 1); i++)
	{
		A[i] /= max1;
		B[i] /= max1; 
	}
	A[n - 1] /= max1;
}


void special_multiply_right (double *a, int n, double * x, int line)
{
	line = line;
	double scal_prod = 0.0;
	for (int i = 0; i < n; i++)
	{
		scal_prod = 0.0;
		for (int j = line; j < n; j++)
		{
			scal_prod += x[j] * a[i * n + j];
		}
		
		for (int j = line; j < n; j++)
		{
			a[i * n + j] += (-2.0 * scal_prod * x[j]);
		}
	}
	
}

void special_multiply_left (double *a, int n, double * x, int line)
{
	line = line;
	double scal_prod = 0.0;
	for (int i = 0; i < n; i++)
	{
		scal_prod = 0.0;
		for (int j = line; j < n; j++)
		{
			scal_prod += x[j] * a[j * n + i];
		}
		for (int j = line; j < n; j++)
		{
			a[j * n + i] += (-2.0 * scal_prod * x[j]);
		}
	}
	
}






void produce_vec(double *a, int n, int line, double *vec)
{
	for (int i = line; i < n; i++)
	{
		vec[i] = a[n * i + (line - 1)];
	}
	
	vec[line] -= norm (vec + line, n - line);
	double N = norm (vec + line, n - line);
	
	if (N > EPS)
	{
		for (int i = line; i < n; i++)
			vec[i] /= N;
	}
}

double norm (double* vec, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++)
    {
        s += vec[i] * vec[i];
    }
    return sqrt (s);
}




int mirror(double *a, int n)
{	
	double * vec = new double [n];
	
	for(int i = 1; i < (n - 1); i++)
	{
		produce_vec (a, n, i, vec);
		special_multiply_left (a, n, vec, i);
		special_multiply_right (a, n, vec, i);
	}		
	
	//print_mat(a, n);

    delete[] vec;
    return 1;
}

void print_mat (double *a, int n)
{	
	int n1=10;
	if(n < 10)
		n1=n;
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n1; j++)
		{
			printf("%f ", a[i * n + j]);
		}
		printf("\n");
	}
}

int read_file (double *a, int n,char ** argv)
{
	FILE *f;
	f = fopen (argv[4], "r");
	if (f == NULL)
	{
		printf("Can't open file\n");
		return -1;
	}
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (fscanf (f, "%lf", &a[i * n + j]) != 1)
			{
				printf ("Can't read from file\n");
				fclose (f);
				return -2;
			}
		}
	}
	fclose (f);
	return 1;
}

double func (int i, int j, int n)
{
	if (i == j && i != n - 1)
	{
		return 1.0;
	}
	else
	{
		if (i == n -1)
		{
			return j + 1;
		}
		
		if (j == n - 1)
		{
			return i + 1;
		}
		return 0.0;
	}
}

void read_func (double *a, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i * n + j] = func (i, j, n);
		}
	}
}

