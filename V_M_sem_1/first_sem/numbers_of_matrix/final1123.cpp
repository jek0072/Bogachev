#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define EPS 1e-15
int read (double *A, char* file, int n);
void func_mat (double *A, int n);
void print_mat (double *A, int n);
void fprint_mat (double *A, int n);
void rotation (double *a, int n);
void determ (double *A, int n, double* cos, double *sin, int i, int j);
void prod_rot_on_mat (double *A, int n, double cos, double sin, int i, int j);
void prod_mat_on_rot (double *A, int n, double cos, double sin, int i, int j);
void residual (double *values, double *a, int n);
void prod_mat (double *a, int a1, int a2, double *b, int b1, int b2, double *c);
double norm (double *A, int n);
void special_RL_prod (double *lr, int n, double *c, int n1);
double* solution (double * A, int n, double eps);
int LR (double *b, double *he, int n, int k);

int LR (double* b, double *he, int n, int k)
{
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
                    he[i * n + j] = b[i * n + j];

    }

    for (int i = 1; i < k; i++)
    {
        //print_mat (b , n);
        if (fabs (he[(i - 1) * n + (i - 1)]) < EPS ){
            //printf ("EROR, can't do LR decomposition\n");
            return -1;
        }

        he[i * n + i - 1] = he[i * n + i - 1] / he[(i - 1) * n + i - 1];
        if (fabs (he[i * n + i - 1]) < 1e-20)
              he[i * n + i - 1] = 0.0;
        for (int j = i; j < k; ++j)
        {
            he[i * n + j] -= he[i * n + i - 1] * he[(i - 1) * n + j];
            if (fabs (he[i * n + j]) < 1e-20)
                he[i * n + j] = 0.0;
        }
    }

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
                    b[i * n + j] = he[i * n + j];

    }
    return 1;
}



void special_RL_prod(double* a, int n, int k)
{


    for (int i = 0; i < k; ++i)
    {
        if (i != 0)
        {
            a[i * n + i - 1] *= a[i * n + i];
            if(fabs (a[i * n + i - 1]) < 1e-20)
                a[i * n + i - 1] = 0.0;
        }
        for (int j = i; j < k - 1; ++j)
        {
            a[i * n + j] += a[i * n + j + 1] * a[(j + 1) * n + j];
            if (fabs (a[i * n + j])< 1e-20)

               a[i * n + j] = 0.0;

        }
    }
}









double* solution (double * A, int n, double eps)
{
        //print_mat(A, n);
    double * own_num = new double [n];
    if (own_num == NULL)
    {
        printf ("ERROR of memory\n");
        return NULL;

    }
    double N = norm (A, n);
    N = N;
    double trace = 0, trace_1 = 0;
    for (int j = 0; j < n; j++)
        trace += A[j * n + j];
    int iterations = 0;
    double invar1 = 0.0,invar2 = 0.0;
    for (int i = 0; i < n; ++i)
    {

        for (int j = 0; j < n; ++j)
            invar1 -= A[i * n + j] * A[j * n + i];
    }
    double te;
    double a, b, c;
    double *B = new double [n * n];
    if (B == NULL)
    {
        printf("ERROR of memory\n");
        delete[] own_num;
        return NULL;

    }

    for (int i = n - 1; i > 1; i--)
    {
        while (fabs(A[(i * n + i - 1)]) > eps)
        {
            //print_mat (A, n);
            //printf ("%e %e\n",fabs(A[(i * n + i - 1)]), eps);
            a = 1;
            b = -(A[(i - 1) * n + (i -1)] + A[i * n + i]);
            c = (A[(i - 1) * n + (i -1)] * A[i * n + i] - A[(i - 1) * n + i] * A[i * n + i - 1]);
            //bool temp;



            if (b * b - 4 * a *c < 0)
            {
                te = fabs(b * b - 4 * a *c);
            }
            else
            {
                double x1 = (-b - sqrt (b * b - 4 * a * c)) / (2 * a);
                double x2 = (-b + sqrt (b * b - 4 * a * c)) / (2 * a);
                if (fabs(A[i * n + i] - x1) <= fabs(A[i * n + i] - x2))
                    te = x1;
                else
                    te = x2;

                //printf("HELLO %d\n", i);
            }
        //te/=100.;
            for (int j = 0; j <= i; j++)
                A[j * n + j] -= te;
            while (LR (A, B, n, i +1) < 0)
            {
                for (int j = 0; j <= i; j++)
                    A[j * n + j] += te;
                te += te / 10.3;
                for (int j = 0; j <= i; j++)
            A[j * n + j] -= te;

                iterations++;
            }

        if (fabs (te) < EPS)
            te = A[i * n + i];



            special_RL_prod (A, n, i + 1);


            for (int j = 0; j <= i; j++)
                A[j * n + j] += te;
            iterations++;



            //printf ("RESIDUAL_1: %e   %d\n", trace_1 - trace, i);
            //printf ("RESIDUAL_2: %e   %d\n", invar2 - invar1, i);
            //print_mat (A, n);

        }
        own_num[i] = A[n * i + i];

    }
    own_num[0] = A[0];

    printf ("\tITERATIONS: %d\n\n", iterations);
    if (n > 1)
    {
        a = 1;
        b = -(A[0] + A[n + 1]);
        c = (A[0] * A[n + 1] - A[1] * A[n]);
    if (b * b - 4 * a * c > 0)
    {
        own_num[0] = (-b - sqrt (b * b - 4 * a * c)) / (2 * a);
        own_num[1] = (-b + sqrt (b * b - 4 * a * c)) / (2 * a);
    }
    else
    {
        if (fabs(b * b - 4 * a * c) < 1e-13)
        {
            own_num[0] = (-b) / (2 * a);
            own_num[1] = (-b) / (2 * a);


        }
        else
        {
            delete[] B;
            delete[] own_num;
            printf ("bad diskr\n");
            printf ("DISKR: %e\n", b*b - 4*a*c);
            printf ("RESIDUAL_1: %e\n", trace_1 - trace);
            invar2 = 0.0;
            for (int i1 = 0; i1 < n; ++i1)
            {

                for (int j = 0; j < n; ++j)
                    invar2 -= A[i1 * n + j] * A[j * n + i1];
            }
            printf ("RESIDUAL_2: %e\n", invar2 - invar1);
           // print_mat (A, n);
            return NULL;
        }
    }
    }
    delete[] B;
    return own_num;
}



double norm (double *A, int n)
{
    double max = 0;
    double t = 0;
    for (int i = 0; i < n; i++)
    {
        t = 0;
        for (int j = 0; j < n; j++)
        {
            t += fabs (A[i * n + j]);
        }
        if (t > max)
            max=t;
    }
    return max;
}



void prod_mat (double *a, int a1, int a2, double *b, int b1, int b2, double *c){
    for (int i = 0; i < a1; i++)
    {
        for (int j = 0; j < b2; j++)
        {
            double s = 0;
            for (int k = 0; k < a2; k++)
            {
                s += a[i * a2 + k] * b[k * b2 + j];
            }
            c[i * b2 + j] = s;
        }
    }
    return;
    return (void)b1;
}

void residual(double *values, double *a, int n)
{
    double inv1 = 0.0;
    double inv2 = 0.0;
    int i, j;
    for (i = 0; i < n; ++i)
    {
        inv1 -= a[i * n + i];
        for (j = 0; j < n; ++j)
            inv2 -= a[i * n + j] * a[j * n + i];
    }

    for (i = 0; i < n; ++i)
    {
        inv1 += values[i];
        inv2 += values[i] * values[i];
    }
    printf("First_RESIDUAL = %e\n", inv1);
        printf("Second_RESIDUAL = %e %d\n", inv2, n);


}
void rotation(double* a, int n)
{


    double x, y, r, cos, sin;

    for (int i = 1; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            x = a[i * n + i - 1];
            y = a[j * n + i - 1];

            if (fabs(y) < 1e-30)
                continue;

            r = sqrt(x * x + y * y);

            if (r < 1e-30)
                continue;

            cos = x / r;
            sin = -y / r;

            a[i * n + i - 1] = r;
            a[j * n + i - 1] = 0.0;

            for (int k = i; k < n; k++)
            {
                x = a[i * n + k];
                y = a[j * n + k];

                a[i * n + k] = x * cos - y * sin;
                a[j * n + k] = x * sin + y * cos;
            }

            for (int k = 0; k < n; k++)
            {
                x = a[k * n + i];
                y = a[k * n + j];

                a[k * n + i] = x * cos - y * sin;
                a[k * n + j] = x * sin + y * cos;
            }
        }
    }
}


int read (double *A, char* file, int n)
{
    FILE *f;
    f = fopen (file, "r");
    if (f == NULL)
    {
        printf ("ERROR, can't open file\n");
        return -1;
    }
    double temp = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (fscanf (f, "%lf", &temp) != 1)
            {
                printf ("ERROR, can't read element\n");
                return -2;
            }
            A[ i * n + j] = temp;
        }
    }

    return 1;
}




void func_mat2 (double *A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[ i * n + j] = 1.0 / (i + j + 1);
        }

    }
}
void func_mat1 (double *A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[ i * n + j] = 0.0;
        }
        A[ i * n + i] = 2;

    }
    for (int i = 0; i < n - 1; i++)
    {
        A[i * n + i + 1] = -1.0;
        A[(i + 1) * n + i] = -1.0;
    }
}
void func_mat	 (double *A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[ i * n + j] = 0.0;
        }


    }

    for (int i = 0; i < n - 1; i++)
    {
        A[i * n + i] = 1.0;
    }
    for (int i = 0; i < n; i++)
    {
        A[i * n + n - 1] = i + 1;
        A[(n - 1) * n + i] = i + 1;
    }

}

void fprint_mat (double *A, int n)
{
    FILE *f;
    f = fopen ("in", "w");
    fprintf (f, "** %d\n", n);
    fprintf (f, "hahahaha\n");

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fprintf (f, "%f ", A[i * n + j]);
        }
        fprintf (f, "\n");
    }
    fprintf (f, "**********************\n\n");
}
void print_mat (double *A, int n)
{
    //printf ( "** %d\n", n);
    //printf ( "hahahaha\n");

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf ("%e ", A[i * n + j]);
        }
        printf ("\n");
    }
    printf ("**********************\n\n");
}


int main (int argc, char **argv){
    if (!(argc == 3 || argc == 4))
    {
        printf ("USAGE: ./a.out  matr_size  Eps\n or\n./a.out  matr_size  Eps  file_name\n");
        return -1;
    }
    int n;
    double eps = 0.0;
    if (sscanf (argv[1], "%d", &n) !=1)
    {
        printf ("USAGE: ./a.out  matr_size  Eps\n or\n./a.out  matr_size  Eps  file_name\n");
        return -2;
    }

    if (sscanf (argv[2], "%lf", &eps) != 1)
    {
        printf ("USAGE: ./a.out  matr_size  Eps\n or\n./a.out  matr_size  Eps  file_name\n");
        return -2;
    }
    double *A = new double [n * n];
    if (A == NULL)
    {
        printf ("ERROR of memory\n");
        return -3;
    }

    if (argc == 4)
    {
        if (read (A, argv[3], n)<0)
        {
            printf ("ERROR of reading\n");
            delete[] A;
            return -4;
        }

    }
    else
    {
        //printf ("here\n");

        func_mat (A, n);
    }

    //double N = norm (A, n);
    //N *= 10;
    for (int i = 0; i < fmin (7, n); i++)
    {
        for (int j = 0; j <fmin (7,n); j++)
            printf ("%f ",A[i * n + j]);
        printf ("\n");
    }

    double invar = 0.0, invar2 = 0.0;
    for (int i = 0; i < n; ++i)
    {

        for (int j = 0; j < n; ++j)
            invar -= A[i * n + j] * A[j * n + i];
    }
    double check = 0;
    for (int i = 0; i < n; i++)
        check += A[i * n + i];

    double t = clock ();
    //printf ("here\n");
    rotation (A, n);
    //print_mat (A, n);
    t = clock () - t;
    double *nums;
    printf ("TIME_ROTATION: %f\n", t / CLOCKS_PER_SEC);
    for (int i = 0; i < n; i++)
        check -= A[i * n + i];
    printf ("RESIDUAL_AFTER_ROTATION_1: %e\n", check);

    for (int i = 0; i < n; ++i)
    {

        for (int j = 0; j < n; ++j)
            invar2 -= A[i * n + j] * A[j * n + i];
    }

    printf ("RESIDUAL_AFTER_ROTATION_2: %e\n", invar2 - invar);

    t = clock();
    nums = solution (A, n, eps);
    if (nums == NULL)
    {
        printf("ERROR in solution\n");
        delete[] A;
        return -13;

    }
    t = clock() - t;
    printf ("TIME_SOLUTION: %f\n", t / CLOCKS_PER_SEC);

    if (argc == 4)
    {
        if (read (A, argv[3], n) < 0)
        {
            printf ("ERROR of reading\n");
            delete[] A;

            delete[] nums;
            return -4;
        }

    }
    else
    {
        //printf ("here\n");

        func_mat (A, n);
    }


    for (int i = 0; i < fmin (n,10); i++)
        printf ("%e\n", nums[i]);
    residual (nums, A, n);
    printf ("****************************************************************************************\n\n\n");
    delete[] nums;
    delete[] A;
    return 0;
}
