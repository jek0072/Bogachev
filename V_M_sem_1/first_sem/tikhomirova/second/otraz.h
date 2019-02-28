#define EPS 1e-15
void read_func (double *a, int n);
double func (int i, int j, int n);
int read_file (double *a, int n,char ** argv);
void print_mat (double *a, int n);
int mirror (double *a, int n);
double norm (double* vec, int n);
void produce_vec(double *a, int n, int line, double *vec);
void special_multiply_left (double *a, int n, double * x, int line);
void special_multiply_right (double *a, int n, double * x, int line);
void change_matr (double *A, double *B, int n);
int n_l (double *A, double *B, int n, double l);
void copy_array (const double *from, double *a,int n);
int solution (double *main, double * sub_main, int n, int k, double eps);

