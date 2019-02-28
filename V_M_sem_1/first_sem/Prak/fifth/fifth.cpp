#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<string.h>
#include<math.h>
# define EPS 1e-15
struct inform
{
    double *a;
    double *column;
    double *string;
    int n1;
    int n2;
    int t;
    int id;
    int *err;
    pthread_mutex_t *mutex;
    pthread_barrier_t *barrier;

};
void* func (void* T)
{
    inform *inf = (inform*)T;
    double *a = inf->a;
    double *column = inf->column;
    double *string = inf->string;
    int n1 = inf->n1;
    int n2 = inf->n2;
    int t = inf->t;
    int id = inf->id;
    string = string;
    a = a;
    n1 = n1;
    n2 = n2;
    column = column;
    t = t;
    id = id;
    //printf ("t= %d\n", t);

    for (int i = id; i < n2; i +=t)
    {
        //printf ("i= %d\n", i);
        double sum = 0.0;
        int num = 0;
        for (int j = 0; j < n1; j++)
        {
            sum += a[j * n2 + i];
            num++;
        }
        column[i] = sum / num;
    }

    for (int i = id; i < n1; i += t)
    {
        double min = a[i * n2 + 0];
        for (int j = 1; j < n2; j++)
        {
            if (min > a[i * n2 + j])
                min = a[i * n2 + j];

        }
        string[i] = min;
    }
    pthread_barrier_wait (inf->barrier);

    for (int i = id; i < n1; i += t)
    {

        for (int j = 0; j < n2; j++)
        {
            if (fabs (a[i * n2 + j] - string[i]) < EPS)
            {
               a[i * n2 + j] = column[j];
            }

        }
    }

    return NULL;
}
int main (int argc, char* argv[])
{
    int n1, n2, t;
    char file_name[100];
    if (argc != 5 )
    {
        printf("USAGE: ./a.out  n1  n2  t  input.txt\n");
        return -1;

    }
    strcpy (file_name, argv[4]);
    if (sscanf (argv[1], "%d", &n1) != 1)
    {
        printf("USAGE: ./a.out  n1  n2  t  input.txt\n");
        return -2;

    }
    if (sscanf (argv[2], "%d", &n2) != 1)
    {
        printf("USAGE: ./a.out  n1  n2  t  input.txt\n");
        return -3;

    }
    if (sscanf (argv[3], "%d", &t) != 1)
    {
        printf("USAGE: ./a.out  n1  n2  t  input.txt\n");
        return -4;

    }
    FILE *f;
    if (n1 <= 0 || n2 <=0 || t <= 0)
    {
        printf("USAGE: ./a.out  n1  n2  t  input.txt\n");
        return -41;

    }
    f = fopen (file_name, "r");
    if ( f == NULL)
    {
        printf ("FILE not exist\n");
        return -5;
    }


    double *ar = new double [n1 * n2];
    if (ar == NULL)
    {
        printf ("ERROR of memory\n");
        return -6;

    }

    for (int i = 0; i < n1 * n2; i++)
    {
        if (fscanf (f, "%lf", &ar[i]) != 1)
        {
            printf("No enough data\n");
            delete[] ar;
            return -7;

        }
    }
    /*for (int i = 0; i < n1; i++)
    {

        printf ("RES: ");
        for (int j = 0; j < n2; j++)
        {
            printf ("%f ", ar[i * n2 + j]);
        }
        printf ("\n");

    }
    printf ("RES:\n");
    */
    if (f != NULL)
        fclose (f);
    printf ("%d %d %d\n",n1,n2,t);
    pthread_t *id = new pthread_t [t];
    if (id == NULL)
    {
        printf ("ERROR id\n");
        delete[] ar;
        return -8;

    }

    inform *inf = new inform [t];
    if (inf == NULL)
    {
        printf ("ERROR inf\n");
        delete[] ar;
        delete[] id;
        return -9;

    }
    pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;
    pthread_barrier_init (&barrier, NULL, t);
    int err = 0;
    double * stolb = new double [n2];
    if (stolb == NULL)
    {
        printf ("ERROR stolb\n");
        delete[] ar;
        delete[] id;
        delete[] inf;
        return -10;

    }
    double *string = new double [n1];
    if (string == NULL)
    {
        printf ("ERROR stolb\n");
        delete[] ar;
        delete[] id;
        delete[] inf;
        delete[] stolb;
        return -11;
    }
    for (int i = 0; i < t; i++)
    {
           inf[i].a = ar;
           inf[i].n1 = n1;
           inf[i].n2 = n2;
           inf[i].id = i;
           inf[i].t = t;
           inf[i].string = string;
           inf[i].err = &err;
           inf[i].barrier = &barrier;
           inf[i].mutex = &mut;
           inf[i].column = stolb;
    }
    int tttt;
    for(int i = 0; i < t; i++)
    {
        tttt = pthread_create (id + i, NULL, &func, inf + i);
        if (tttt < 0)
        {
            printf ("ERROR of create\n");
            pthread_barrier_destroy(&barrier);
            delete[] ar;
            delete[] id;
            delete[] inf;
            delete[] stolb;
            delete[] string;
            return -666;
        }
    }
    for (int i = 0; i < t; i++)
        pthread_join (id[i], NULL);

    printf ("sred:\n");
    for (int i = 0; i < n2; i++)
        printf ("%f ", stolb[i]);
    printf ("\n");
    printf ("min:\n");
    for (int i = 0; i < n1; i++)
        printf ("%f ", string[i]);
    printf ("\n");
    for (int i = 0; i < n1; i++)
    {

        printf ("RES: ");
        for (int j = 0; j < n2; j++)
        {
            printf ("%f ", ar[i * n2 + j]);
        }
        printf ("\n");

    }
    printf ("RES:\n");



    pthread_barrier_destroy (&barrier);
    delete[] ar;
    delete[] id;
    delete[] inf;
    delete[] string;
    delete[] stolb;
    return 0;
}
