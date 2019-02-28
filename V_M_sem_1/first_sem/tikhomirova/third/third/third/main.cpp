#include<sys/time.h>
#include<stdio.h>
#include<pthread.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"otraz.h"




int main (int argc, char** argv)
{
	if (!(argc == 3 || argc == 4))
	{
		printf ("USAGE: ./a.out	 size_matr  number_of_threads  file_name\n or\n./a.out  size_matr  number_of_threads\n");
		return -1;
	}
	int n, t;
	
	double *a;
	n = atoi (argv[1]);
	
	t = atoi (argv[2]);
	if (n <= 0 || t <= 0){
		printf("Incorect data\n");
		return -11;
	}
		
	a = new double [n* n];
	
	if (a == NULL)
	{
		printf("ERROR of memory\n");
		return -2;
	}
	
	double *b = new double [n];
	if (b == NULL)
	{
		printf("ERROR of memory\n");
		delete[] a;
		
		return -13;
	}	
	
	if (argc == 3)
	{
		read_func(a, b, n);
	}
	else
	{
		if (read_file(a, b, n, argv) < 0)
		{
			printf("ERROR in the reading function\n");
			delete[] a;
			delete[] b;
			return -3;
		}
			
    }

	pthread_t *id = new pthread_t [t];
	if (id == NULL)
	{
		printf ("ERROR of memory\n");
		delete[] a;
		delete[] b;
		return -14;
	}
		
	inf_thr *INF = new inf_thr [t];
	if (INF == NULL)
	{
		printf ("ERROR of memory\n");
		delete[] a;
		delete[] b;
		delete[] id;
		return -15;
	}
	
	double *shared_vec = new double [n];
	if (shared_vec == NULL)
	{
		printf ("ERROR of memory\n");
		delete[] a;
		delete[] b;
		delete[] id;
		delete[] INF;
		return -16;	
	}
	
	double time = 0;
	pthread_barrier_t barrier;
	pthread_barrier_init (&barrier, NULL, t);
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;	
	
	for (int i = 0; i < t; i++)
	{
		INF[i].a = a;
		INF[i].b = b;
		INF[i].id = i;
		INF[i].n = n;
		INF[i].t = t;
		INF[i].time = &time;
		INF[i].shared_vec = shared_vec;
		INF[i].barrier = &barrier;
		INF[i].mutex = &mutex;
	}
	
	
	int ttt;
    double t0 = get_full_time();
	for (int i = 0; i < t; i++)
	{
		ttt = pthread_create (id + i, NULL, thread_function, INF + i);
        if (ttt < 0)
		{
			printf ("YOU have big problems\n");
			return -666;
		}
		
	}
	
	
	for (int i = 0; i < t; i++)
	{
		pthread_join (id[i], NULL);
	}
    printf ("TIME: %f\n", get_full_time () - t0);
	
	
	print_answer(b, n, 10);
		
	double* x = new double [n];
	for (int i = 0; i < n; i++)
		x[i] = b[i];
	
	if (argc == 3)
	{
		read_func(a, b, n);
	}
	else
	{
		if (read_file(a, b, n, argv) < 0)
		{
			printf("ERROR in the reading function\n");
            pthread_barrier_destroy (&barrier);
            delete[] x;
            delete[] a;
            delete[] b;
            delete[] id;
            delete[] INF;
            delete[] shared_vec;
			return -3;
		}
			
	}

	printf ("residual: %e\n", residual (a, b, x, n));	
	 
	
	
			
	pthread_barrier_destroy (&barrier); 	
	delete[] x;
	delete[] a;
	delete[] b;
	delete[] id;
	delete[] INF;
	delete[] shared_vec;
	return 0;
}
