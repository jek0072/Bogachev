#include<stdio.h>
#include<math.h>
#include<sys/time.h>
#include<stdio.h>
#include<pthread.h>

#include<stdlib.h>
#include<time.h>
#include"otraz.h"


double get_full_time(){
	struct timeval buf;
	gettimeofday(&buf,0);
	return (double)buf.tv_sec+(double)buf.tv_usec/1000000.0;
}

double residual (double *a, double *b, double *x, int n){
    double *res=new double[n];    
    for(int i=0;i<n;i++){
        double s=0.0;
        for(int j=0;j<n;j++){
            s+=a[i*n+j]*x[j];
        }
        res[i]=s;
        res[i]-=b[i];
    }
    double p=0.0;
    for(int i=0;i<n;i++)
        p+=res[i]*res[i];
    delete[] res;
    return sqrt(p);
}






void special_multiply_vec (double *b, int n, double * x, int line)
{
	line = line;
	double scal_prod = 0.0;
	
	scal_prod = 0.0;
	
	for (int j = line; j < n; j++)
	{
		scal_prod += x[j] * b[j];
	}
	
	for (int j = line; j < n; j++)
	{
		b[j] += (-2.0 * scal_prod * x[j]);
	}
	
	
}

void* thread_function (void *T)
{
	inf_thr *inf = (inf_thr*)T;
	double *a = inf -> a;
	double *b = inf -> b;
	pthread_mutex_t *mutex = inf -> mutex;
	pthread_barrier_t *barrier = inf -> barrier;
	int n = inf -> n;
	int id = inf -> id;
	int t = inf -> t;
	double *shared_vec = inf -> shared_vec;
	double *time = inf -> time;
	double * vec;
    vec = new double [1];
    delete[] vec;
	if (id == 0)
    {
        vec = new double [n];

    }

	for(int i = 0; i < (n - 1); i++)
	{
		if (id == 0)
		{
			produce_vec (a, n, i, vec);
			for (int j = i; j < n; j++)
				shared_vec[j] = vec[j];
		}
		pthread_barrier_wait (barrier);
		
		special_multiply_left (a, n, shared_vec, i, id, t);
		
		if(id == 0)
			special_multiply_vec (b, n, vec, i);
		
		pthread_barrier_wait (barrier);
	}
		
	for (int i = 0; i < n; i++){
        if (fabs (a[(n - 1 - i) * n + (n - 1 - i)]) < EPS)
        {	
        	pthread_mutex_lock (mutex);
            *time = -1;
        	pthread_mutex_unlock (mutex);

            pthread_barrier_wait(barrier);

            return NULL;
        }
        double tt = b[n - 1 - i] / a[(n - 1 - i) * n + (n - 1 - i)];
        a[(n - 1 - i) * n + (n - 1 - i)] = 1.0;
        b[n - 1 - i] = tt;
    
        for (int j = id; j < n - i - 1; j += t)
        {
        	
	        b[j] -= tt * a[j * n + (n - i - 1)];
        }
        pthread_barrier_wait (barrier);
    }	
	if (id == 0)
	{

        delete[] vec;
	}
	return NULL;
}




void special_multiply_left (double *a, int n, double * x, int line, int id, int t)
{
	line = line;
	double scal_prod = 0.0;
	for (int i = n - 1 - id; i >= line; i -= t)
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
		vec[i] = a[n * i + (line)];
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

int read_file (double *a, double *b, int n, char ** argv)
{
	FILE *f;
	f = fopen (argv[3], "r");
	printf ("%s\n", argv[3]);
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
		
		if (fscanf (f, "%lf", &b[i]) != 1)
		{
			printf ("Error of reading \n");
			fclose(f);
			return -3;
		}
	}
	fclose (f);
	return 1;
}

double func (int i, int j, int n)
{
	n = n;
	return fabs (i - j);
}

void read_func (double *a, double *b, int n)
{	
	for (int i = 0; i < n; i++)
	{
		b[i] = 0.0;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i * n + j] = func (i, j, n);
			
			if (j % 2 == 1)
				b[i] += a[i * n + j];
		}
	}
}

