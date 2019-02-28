#define EPS 1e-10
void read_func (double *a, double *b, int n);
double func (int i, int j, int n);
int read_file (double *a, double *b, int n,char ** argv);

void print_mat (double *a, int n);

double norm (double* vec, int n);
void produce_vec(double *a, int n, int line, double *vec);
void special_multiply_left (double *a, int n, double * x, int line, int id, int t);
int inverse_gauss_for_vector (double *A, double *vec, int n);
void print_answer(double* vec,int n,int m);
void special_multiply_vec (double *b, int n, double * x, int line);
double residual (double *a, double *b, double *x, int n);
double get_full_time();
void* thread_function (void *T);


struct inf_thr
{
	double *a;
	double *b;
	pthread_mutex_t *mutex;
	pthread_barrier_t *barrier;
	int n;
	int id;
	int t;
	double *shared_vec;
	double *time;
};


