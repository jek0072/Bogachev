#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<string.h>
#include<math.h>
# define EPS 1e-15
struct inform
{
    double *a;
    int n1;
    int n2;
    int t;
    int id;
    int *err;
    pthread_mutex_t *mutex;
    pthread_barrier_t *barrier;
    double *sum;
    int *amount;
    int beg;
    int end;

};
void* func (void* T)
{
    inform *inf = (inform*)T;
    double *a = inf->a;
    int n1 = inf->n1;
    int n2 = inf->n2;
    int t = inf->t;
    int id = inf->id;
    int beg = inf->beg;
    int end = inf->end;
    a = a;
    n1 = n1;
    n2 = n2;
    t = t;
    id = id;
    beg = beg;
    end = end;
    double sum = 0.0;
    int am = 0;
    for (int i = beg; i < end; i++)
    {
        for (int j = 0; j < n2; j++)
        {
            if (i - 1 >= 0 && i + 1 < n1 && j - 1 >= 0 && j + 1 < n2)//check 1
            {
                if (a[i * n2 + j ] > a[i * n2  + j + 1] && a[i * n2 + j ] > a[i * n2  + j - 1] && a[i * n2 + j ] > a[(i + 1) * n2  + j] && a[i * n2 + j ] > a[(i - 1) * n2  + j])// check 2
                {
                    sum += a[i * n2 + j];
                    am++;
                }
            }
        }

    }
    printf ("ID: %d\t%d\t%d\n", id, beg, end);
    pthread_barrier_wait (inf->barrier);
    pthread_mutex_lock (inf->mutex);
    *inf->amount += am;
    *inf->sum += sum;
    pthread_mutex_unlock (inf->mutex);

    pthread_barrier_wait (inf->barrier);
    if (*inf->amount == 0)
    {
        pthread_mutex_lock (inf->mutex);
        *inf->err = -1;
        pthread_mutex_unlock (inf->mutex);
        return NULL;
    }
    if (id == 0)
        *inf->sum /= *inf->amount;
    pthread_barrier_wait (inf->barrier);

    int num_string_in_thread = end - beg;
    double *Beg = new double [n2];
    double *End = new double [n2];
    Beg = Beg;
    End = End;
    if (num_string_in_thread > 2)
        num_string_in_thread = 2;

    if (num_string_in_thread == 1)
    {
        if (beg - 1 >=0 && beg + 1 < n1)
        {
            for (int i = 0; i < n2; i++)
            {
                if (i - 1 >= 0 && i + 1 < n2)
                {
                    if(a[beg * n2 + i ] > a[beg * n2  + i + 1] && a[beg * n2 + i ] > a[beg * n2  + i - 1] && a[beg * n2 + i ] > a[(beg + 1) * n2  + i] && a[beg * n2 + i ] > a[(beg - 1) * n2  + i])
                    {
                        Beg[i] = *inf->sum;
                    }
                    else
                        Beg[i] = a[beg * n2 + i];
                }
                else
                    Beg[i] = a[beg * n2 + i];
            }

        }
        else
        {
            for (int i = 0; i < n2; i++)
            {
                Beg[i] = a[beg * n2 + i];
            }
        }

    }
    if (num_string_in_thread == 2)
    {
        if (beg - 1 >=0 && beg + 1 < n1)
        {
            for (int i = 0; i < n2; i++)
            {
                if (i - 1 >= 0 && i + 1 < n2)
                {
                    if(a[beg * n2 + i ] > a[beg * n2  + i + 1] && a[beg * n2 + i ] > a[beg * n2  + i - 1] && a[beg * n2 + i ] > a[(beg + 1) * n2  + i] && a[beg * n2 + i ] > a[(beg - 1) * n2  + i])
                    {
                        Beg[i] = *inf->sum;
                    }
                    else
                        Beg[i] = a[beg * n2 + i];
                }
                else
                    Beg[i] = a[beg * n2 + i];
            }

        }
        else
        {
            for (int i = 0; i < n2; i++)
            {
                Beg[i] = a[beg * n2 + i];
            }
        }




        if ((end - 1) - 1 >=0 && (end - 1) + 1 < n1)
        {
            for (int i = 0; i < n2; i++)
            {
                if (i - 1 >= 0 && i + 1 < n2)
                {
                    if(a[(end - 1) * n2 + i ] > a[(end - 1) * n2  + i + 1] && a[(end - 1) * n2 + i ] > a[(end - 1) * n2  + i - 1] && a[(end - 1) * n2 + i ] > a[((end - 1) + 1) * n2  + i] && a[(end - 1) * n2 + i ] > a[((end - 1) - 1) * n2  + i])
                    {
                        End[i] = *inf->sum;
                    }
                    else
                        End[i] = a[(end - 1) * n2 + i];
                }
                else
                    End[i] = a[(end - 1) * n2 + i];
            }

        }
        else
        {
            for (int i = 0; i < n2; i++)
            {
                End[i] = a[(end - 1) * n2 + i];
            }
        }

    }
    pthread_barrier_wait (inf->barrier);
    double *fir = new double [n2];
    double *sec = new double [n2];
    bool tt = 0;
    for (int i = beg + 1; i < end - 1; i++)
    {

        //printf ("HAHAHAHA %d\n", i);
        for (int j = 0; j < n2; j++)
        {


            if (j - 1 >= 0 && j + 1 < n2 && a[i * n2 + j ] > a[i * n2  + j + 1] && a[i * n2 + j ] > a[i * n2  + j - 1] && a[i * n2 + j ] > a[(i + 1) * n2  + j] && a[i * n2 + j ] > a[(i - 1) * n2  + j])// check 2
            {
                 //printf ("%d %d\n", i, j);
                sec[j] = *inf->sum;
            }
            else
                sec[j] = a[i * n2 + j];


        }
        if (tt != 0)
        {
            //printf ("HAHAHAHA %d\n", i);

            for (int j = 0; j < n2; j++)
            {
                a[(i - 1) * n2 + j] = fir[j];
            }
        }

        for (int j = 0; j < n2; j++)
        {
            fir[j] = sec[j];
        }
        tt = 1;

    }

    if  (tt == 1)
    {
        //printf ("TITITITI\n");
        for (int j = 0; j < n2; j++)
        {
            a[(end - 1 - 1) * n2 + j] = fir[j];
        }

    }
    pthread_barrier_wait (inf->barrier);

    if (num_string_in_thread == 1)
    {    for (int i = 0; i < n2; i++)
        {
            a[beg * n2 + i] = Beg[i];
        }
    }

    if (num_string_in_thread == 2)
    {    for (int i = 0; i < n2; i++)
        {
            a[beg * n2 + i] = Beg[i];
            a[(end - 1) * n2 + i] = End[i];
        }
    }

    delete[]fir;
    delete[] sec;
    delete[] Beg;
    delete[] End;
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

    if (t > n1)
        t = n1;

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
    int amount = 0;
    double sum = 0.0;
    for (int i = 0; i < t - 1; i++)
    {
           inf[i].a = ar;
           inf[i].n1 = n1;
           inf[i].n2 = n2;
           inf[i].id = i;
           inf[i].t = t;
           inf[i].sum = &sum;
           inf[i].amount = &amount;
           inf[i].err = &err;
           inf[i].barrier = &barrier;
           inf[i].mutex = &mut;
           inf[i].beg = i * (n1 / t);
           inf[i].end = (i + 1) * (n1 / t);
    }
    inf[t - 1].a = ar;
    inf[t - 1].n1 = n1;
    inf[t - 1].n2 = n2;
    inf[t - 1].id = t - 1;
    inf[t - 1].t = t;
    inf[t - 1].sum = &sum;
    inf[t - 1].amount = &amount;
    inf[t - 1].err = &err;
    inf[t - 1].barrier = &barrier;
    inf[t - 1].mutex = &mut;
    inf[t - 1].beg = (t - 1) * (n1 / t);
    inf[t - 1].end = n1;
    int tttt;
    /*
    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < n2; j++)
        {
            printf ("%f ", ar[i * n2 +j]);

        }
        printf ("\n");
    }
    */
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
            return -666;
        }
    }
    for (int i = 0; i < t; i++)
        pthread_join (id[i], NULL);
    printf ("res: %f\n", sum);
    if( err == -1)
    {
        printf ("HELLO\n");
    }
    else
    {
        printf ("RES:\n");
        for (int i = 0; i < n1; i++)
        {

            printf ("RES: ");
            for (int j = 0; j < n2; j++)
            {
                printf ("%f ", ar[i * n2 + j]);

            }
            printf ("\n");        }

        printf ("RES:\n");
        printf ("hehe\n");

    }


    pthread_barrier_destroy (&barrier);
    delete[] ar;
    delete[] id;
    delete[] inf;
    return 0;
}
