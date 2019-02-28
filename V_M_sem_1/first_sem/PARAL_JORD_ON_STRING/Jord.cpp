#include<pthread.h>
#include <sched.h>	
#include <sys/resource.h>
#include<stdio.h>
#include<stdlib.h>
#include<ctime>
#include<sys/time.h>
#include<cmath>
#include <unistd.h>

#define EPS 1e-10
//in.txt n d     or      n d
//-ffast-math -O3 -march=native
int read(double*mat,FILE *f,int n,int d);//OK
void print_mat(double *a,int n);//OK
void print_convert_mat(double *a,int d,int k,int p);//OK
double* prod_ab(double* a,int a1,int a2, double *b,int b1,int b2,double* c);//OK
double norm(double *a,int n);//OK
int inverse_mat(double* A,double* E,int m,double N);//OK inverse matrix in E      //1-Ok, 0-Err
void a_minus_b(double *a,double *b,int n,int m);//OK
int Jordan(double *A, double *E, int k, int p, int d, int *index_column, int t, pthread_t *id,int argc,char** argv);

int search_block(double *A,double *E,int line,int k,int d,int p,int t,double N,int id);
void copy(double *a,double *b,int n,int m);//OK
void null(double *a,int n);//OK
void unit(double *a,int n);//OK
void convert(double *a,double *b,int d,int p,int k);//convertation from  a in b
void func_mat(double *A,int k, int d, int p);//GILBERT
void PRINT_INVERSE(double *a,int d,int k,int p,int *ind);//PRINT AFTER TRANSFORMATION
int read(double* mat, FILE *f,int n);
void convert_E(double *a,double *b,int d,int p,int k,int *ind);
void CONVERTATION_STRING(int *ind,int n);
void* thread_0(void* T);
int typical_prod(double *a,int a1,int a2, double *b,int b1,int b2,double *c);
double residual(double* a, int n, int t);
void* count_residual(void*T);


void put(double*A,double *res,int num_st,int k,int n);
void get_string(double*A,double *res,int num_st,int k,int n);

void get_square(double *A,double *res,int num_st,int num_col,int k,int n);
double get_full_time();
double get_time();

double get_full_time(){
	struct timeval buf;
	gettimeofday(&buf,0);
	return (double)buf.tv_sec+(double)buf.tv_usec/1000000.0;
}	
double get_time(){
	struct rusage buf;
	getrusage(RUSAGE_THREAD,&buf);
	return (double)buf.ru_utime.tv_sec+(double)buf.ru_utime.tv_usec/1000000.;
}


struct inform_thread{
    double *A;
    double *E;
    int* ind_col;
    double *N;
    double *TIME;
    int p;
    int k;
    int t;
    int d;
	int argc;
	char** argv;
    int *column;
    double* inverse_share;
    double* norm;
    pthread_mutex_t *inverse_mutex;
	double *res;
    int* err;
    pthread_mutex_t *error;

    pthread_barrier_t *barrier;
    int id;
    pthread_mutex_t* help_mutex;
};

struct inf_res{
    double *A;
    double *Inv;
    int n;
    int id;
    int t;
    double *res;
    pthread_barrier_t *barrier;
};

void get_square(double *A,double *res,int num_st,int num_col,int k,int n){
	for(int i=0;i<k;i++){
		for(int j=0;j<k;j++){
			res[i*k+j]=A[n*(i+k*num_st)+(j+k*num_col)];
		}
	}
}

void get_string(double*A,double *res,int num_st,int k,int n){
    for(int i=0;i<k;i++){
        for(int j=0;j<n;j++){
            res[i*n+j]=A[n*(i+num_st*k)+(j)];
        }
    }
}

void put(double*A,double *res,int num_st,int k,int n){//block string
    for(int i=0;i<k;i++){
        for(int j=0;j<n;j++){
            A[n*(i+num_st*k)+(j)]+=res[i*n+j];
        }
    }
}

void* count_residual(void*T){
    inf_res *inf=(inf_res*)T;
    int k=inf->n/inf->t;
    int n=inf->n;
    int t=inf->t;
    int id=inf->id;
    double* A=inf->A;
    double *Inv=inf->Inv;
    double* res=inf->res;
    int j=id;
    
    double *c=new double[n*n];
    double *b=new double[k*n];
    double *a=new double[k*k];
    	

    double t0=get_full_time();
    for(int i=0;i<t;i++){
        if(j>=t)
            j=0;
        get_square(A,a,id,j,k,n);
        get_string(Inv,b,j,k,n);
        prod_ab(a,k,k,b,k,n,c);
        put(res,c,id,k,n);
        //pthread_barrier_wait(inf->barrier);
        /*if(id==0){
		    printf("\nHERE\n");
		    print_mat(res,n);
		    printf("\n\n");
        }*/
        j++;
    }
    pthread_barrier_wait(inf->barrier);
        
    if(id == 0 && (n % t) != 0)
    {
    
    	for(int i = 0; i < n; i++)
    	{
    		for(int j = 0; j < n; j++)
    		{
    			double sum=0.0;
    			
				for(int k1 = k * t; k1 < n; k1++)
				{
					sum += A[i * n + k1] * Inv[k1 * n + j];	
				}
				
				res[i * n + j] += sum;
				
			}	
    	}
    	
    	prod_ab( A + n * (k * t), n - k * t, n - (n % t), Inv, n - (n % t), n, c);
    	
    	for(int i = k * t; i < n; i++)
    	{
    		//printf("here\n");
    		for(int j = 0; j < n - (n % t); j++)
    		{
    			res[i * n + j] += c[(i - k* t) * n + j];
    		}
    	}
    		
    }
    
    t0=get_full_time()-t0;
    if (id == 0) printf (" Elapsed_RESIDUAL = %.5f\n", t0);
    
    //if(id==0)
    	//print_mat(res,n);
    delete[] a;
	delete[] b;
	delete[] c;
    return NULL;


}

double residual(double* a,double *inv, int n, int t){
    pthread_t *id=new pthread_t [t];
    inf_res *inf=new inf_res[t];
    double *res=new double[n*n];
    for(int i=0;i<n;i++){
    	for(int j=0;j<n;j++){
    		res[i*n+j]=0.0;
    	}
    }
    pthread_barrier_t b;
    pthread_barrier_init(&b,NULL,t);
    for(int i=0;i<t;i++){
        inf[i].A=a;
        inf[i].Inv=inv;
        inf[i].t=t;
        inf[i].id=i;
        inf[i].res=res;
        inf[i].n=n;
        inf[i].barrier=&b;
    }
    int t_t;
    for(int i=0;i<t;i++){
        t_t=pthread_create(id+i,NULL,&count_residual,inf+i);
        if(t_t<0){
        	printf("ahahah\n");
        	
			delete[]id;
			delete[] inf;
			delete[]res;
			pthread_barrier_destroy(&b);
        	return -1;
        }
    }
    for(int i=0;i<t;i++){
        pthread_join(id[i],NULL);

    }
    int k=n/t;
    for(int i=t*k;i<n;i++){//peace we need to prod to the string
    	for(int j=0;j<n;j++){
    		double s=0.0;
    		for(int k1=0;k1<n;k1++){
    			s+=a[i*n+k1]*inv[k1*n+j];
    		}
    		res[i*n+j]=s;
    	}
    }
    
    for(int i=0;i<t*k;i++){
    	for(int j=t*k;j<n;j++){
    		double s=0.0;
    		for(int k1=0;k1<n;k1++){
    			s+=a[i*n+k1]*inv[k1*n+j];
    		}
    		res[i*n+j]=s;
    	}
    	
    }
    //printf("\n");
    //print_mat(res,n);
    
    
    for(int i=0;i<n;i++){
        res[i*n+i]-=1;
    }
    //
    double r=norm(res,n);
    
    delete[]id;
	delete[] inf;
	delete[]res;
    pthread_barrier_destroy(&b);
    printf("RESIDUAL: %e\n",r);
    return r;

}

void* thread_0 (void* T){
    inform_thread *inform = (inform_thread*)(T);
    int d = inform -> d;
    int p = inform -> p;
    int k = inform -> k;
    //int n=k*d+p;
    int argc = inform -> argc;
    char** argv = inform -> argv;
    int t = inform -> t;
    int id = inform -> id;
    FILE* f;
    double* Ed = new double [d*d];
    double* TEMP = new double [d*d];
    double* E = inform -> E;
    double *A = inform -> A;
    int res = 0;
    //printf("id=%d \n",id);
		
	//printf("HELLO\n");
	//printf("d= %d  k=%d p=%d t=%d\n",d,k,p,t);
	
	for (int i = k - 1 - id; i >= 0; i -= t){
		for (int j = 0; j < k; j++){
			for (int i1 = 0; i1 < d; i1++){
				for (int i2 = 0; i2 < d; i2++){
					A[i*(d*d*k+p*d)+j*d*d+i1*d+i2]=0.0;
					E[i*(d*d*k+p*d)+j*d*d+i1*d+i2]=0.0;
				
				}
			}
		}
	}
	
	for (int i = k - 1 - id; i >= 0; i -= t){
		for (int i1 = 0; i1 < d; i1++){
			for (int i2 = 0; i2 < p; i2++){
				A[i*(d*d*k+p*d)+k*d*d+i1*p+i2]=0.0;
				E[i*(d*d*k+p*d)+k*d*d+i1*p+i2]=0.0;
			}
		}
	}
	
	if(id==0){
		for (int j = 0; j < k; j++){
			for (int i1 = 0; i1 < p; i1++){
				for (int i2 = 0; i2 < d; i2++){
					A[k*(d*d*k+p*d)+j*p*d+i1*d+i2]=0.0;
					E[k*(d*d*k+p*d)+j*p*d+i1*d+i2]=0.0;
				}
			}
		}

		
		for (int i1 = 0; i1 < p; i1++){
			for(int i2=0;i2<p;i2++){
				A[k*(d*d*k+p*d)+k*p*d+i1*p+i2]=0.0;
				E[k*(d*d*k+p*d)+k*p*d+i1*p+i2]=0.0;
			}
		}
		
	
	}
	pthread_barrier_wait(inform->barrier);

	if(id==0){
		//print_convert_mat(A,d,k,p);
		if(argc==5){
		    f = fopen(argv[4], "r");
			printf("%s\n",argv[4]);
					
		    if(f == NULL){                                             //bad file
		        printf("CAN'T OPEN FILE\n");
				
		        
		        pthread_mutex_lock(inform->error);
		        	*inform->err=-2;
		        pthread_mutex_unlock(inform->error);
		        
		        delete[] TEMP;
		        delete[] Ed;
			    pthread_barrier_wait(inform->barrier);
		        return NULL;
		    }

		    if(f!=NULL && (read(A,f,k*d+p,d))<0){
		        printf("ERROR of reading or \n");
		        if(f!=NULL)
		        	fclose(f);
		        
		        pthread_mutex_lock(inform->error);
		        	*inform->err=-2;
		        pthread_mutex_unlock(inform->error);
		        delete[] TEMP;
		        delete[] Ed;
			    pthread_barrier_wait(inform->barrier);
			    return NULL;
		    }
		    if(f!=NULL)
		    	fclose(f);
		   // print_convert_mat(A,d,k,p);
		}
		else{
		    func_mat(A,k,d,p);
		}
		*inform->N = norm(A,k*d+p);
		//null(Ed,d*k+p);
		for(int i=0;i<k;i++){
		    for(int j1=0;j1<d;j1++){

		        E[i*(d*d*k+p*d)+i*d*d+j1*d+j1]=1.0;

		    }
		}
		
		for(int i=0;i<p;i++){
		    E[k*(d*d*k+p*d)+d*p*k+i*p+i]=1.0;
		}
		//print_convert_mat(A,d,k,p);
		
	}
	    
	pthread_barrier_wait(inform->barrier);
	if(*inform->err==-2){
		printf("HERE\n");
		delete[] TEMP;
        delete[] Ed;
		return NULL;
	}
	//pthread_barrier_wait(inform->barrier);



    double t0 = 0.0, t1 = 0.0;
    t0 = get_time ();
    t1 = get_full_time ();
//		printf("%f\n",*inform->N);
    for (int i = 0; i < k; i++){

        res = search_block (A, Ed, i, k, d, p, t, *inform -> N, id);//returns -1 if cannot find inverse
       
       	
       	
       if(res >= 0)
       {
        	
        
            pthread_mutex_lock(inform->inverse_mutex);
            if (fabs (*(inform -> norm)) < EPS)
            {
                *(inform -> column) = res;
                //printf("HAHAH\n\n\n");
                *(inform -> norm) = norm(Ed,d);
                //printf("HAHAH %lf\n\n\n",*(inform->norm));
            }
            else
            {
                double nn = norm(Ed,d);
                if(*(inform -> norm) > nn)
                {
                    *(inform -> column) = res;
                   // printf("HAHAH %lf\n\n\n",nn);
                    *(inform -> norm) = nn;
                }
            }


            pthread_mutex_unlock(inform->inverse_mutex);

        }
        
        pthread_barrier_wait (inform -> barrier);//WAIT INVERSE
        
        if (id == 0)
        {
		    
		    if (!(fabs(*inform -> norm) <  EPS) )
		    {
		        *inform -> norm = 0.0;//i and res
		        int ttt = inform -> ind_col[i];
		        inform -> ind_col[i] = inform -> ind_col [*inform -> column];
		        inform -> ind_col[*inform -> column] = ttt;
		    }
		    else{
		        *inform->err = -1;
		        *(inform->column) = 1;
		    }
        }
        
		
       //when we multiply string how do this
        for(int j = k - 1 - id; j >= 0;j -= t)
        {
            copy (TEMP, A + j * (d * d * k + p * d) + *inform -> column * d * d, d, d);
            copy (A + j * (d * d * k + p * d) + *inform -> column * d * d, A + j * (d * d * k + p * d) + i * d * d, d, d);
            copy (A + j * (d * d * k + p * d) + i * d * d, TEMP, d, d);
        }
        
        if (p != 0 && inform -> id == 0)
        {
		    copy (TEMP, A + k * (d * d * k + p * d) + *inform -> column * d * p, p, d);
		    copy (A + k * (d * d * k + p * d) + *inform -> column * p * d, A + k * (d * d * k + p * d) + i * p * d, p, d);
		    copy (A + k * (d * d * k + p * d) + i * p * d, TEMP, p, d);
		}
		
        pthread_barrier_wait (inform -> barrier);//WAIT for copy
        
       	if (*inform -> err == -1)
       	{
            delete[] TEMP;
            delete[] Ed;
            return NULL;
        }
        copy (TEMP, A + i * (d * d * k + p * d) + i * d * d, d, d);
        inverse_mat(TEMP,Ed,d,*inform->N);
         
        for (int j = k - 1 - id ; j > i; j -= t)
        {                                         //prod line to inverse element
            prod_ab(Ed,d,d,A+i*(d*d*k+p*d)+j*d*d,d,d,TEMP);
            copy(A+i*(d*d*k+p*d)+j*d*d,TEMP,d,d);
        }


        for (int j = k - 1 - id; j >= 0; j -= t)
        {                                         //prod E line to inverse element
            prod_ab(Ed,d,d,E+i*(d*d*k+p*d)+j*d*d,d,d,TEMP);
            copy(E+i*(d*d*k+p*d)+j*d*d,TEMP,d,d);
        }

        if(p!=0 && id==0){
            prod_ab(Ed,d,d,A+i*(d*d*k+p*d)+k*d*d,d,p,TEMP);
            copy(A+i*(d*d*k+p*d)+k*d*d,TEMP,d,p);

            prod_ab(Ed,d,d,E+i*(d*d*k+p*d)+k*d*d,d,p,TEMP);
            copy(E+i*(d*d*k+p*d)+k*d*d,TEMP,d,p);

        }
        pthread_barrier_wait (inform -> barrier);
        
        
        

        
        
        
        for(int j = k - 1 - id; j >= 0; j -= t)
        {// high correction
            for(int y = k - 1; y > i; y--)
            {
            	if (j != i)
                	a_minus_b(A+j*(d*d*k+p*d)+d*d*y,prod_ab(A+i*d*d+j*(d*d*k+d*p),d,d,A+i*(d*d*k+p*d)+d*d*y,d,d,TEMP),d,d);
        
            }
        }////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for(int j = k - 1 - id; j >= 0; j -= t)
        {// high correction
            for(int y = k - 1; y >= 0; y--)
            {
            	if (j != i)
                	a_minus_b(E+j*(d*d*k+p*d)+d*d*y,prod_ab(A+i*d*d+j*(d*d*k+d*p),d,d,E+i*(d*d*k+p*d)+d*d*y,d,d,TEMP),d,d);
        
            }
        }

		//will it be faster?
        if (p != 0)
        {
        	///////
        	if (id == 0)
        	{
		        for(int j = k - 1; j > i; j--)
		        {//lowest correction
		                a_minus_b(A+k*(d*d*k+p*d)+d*p*j,prod_ab(A+i*p*d+k*(d*d*k+d*p),p,d,A+i*(d*d*k+p*d)+d*d*j,d,d,TEMP),p,d);
		        }
		        a_minus_b(A+k*(d*d*k+p*d)+d*p*k,prod_ab(A+i*p*d+k*(d*d*k+d*p),p,d,A+i*(d*d*k+p*d)+d*d*k,d,p,TEMP),p,p); 
		        
		        for(int j = k - 1; j >= 0; j--)
		        {//lowest correction
		                a_minus_b(E+k*(d*d*k+p*d)+d*p*j,prod_ab(A+i*p*d+k*(d*d*k+d*p),p,d,E+i*(d*d*k+p*d)+d*d*j,d,d,TEMP),p,d);
		        }
		        a_minus_b(E+k*(d*d*k+p*d)+d*p*k,prod_ab(A+i*p*d+k*(d*d*k+d*p),p,d,E+i*(d*d*k+p*d)+d*d*k,d,p,TEMP),p,p);
            
            }
            /////
            
				for(int j = k - 1 - id; j >= 0; j -= t){//right correction
		                if (j != i)
		                	a_minus_b ( A + j * (d * d * k + p * d) + d * d * k, prod_ab (A + i * d * d + j * (d * d * k + d * p), d, d, A + i * (d * d * k + p * d) + d * d * k, d, p, TEMP), d, p);
		        }
            
		    
		    	
            
            for(int j = k - 1 - id; j >= 0; j -= t){//right correction
		                if (j != i)
		                	a_minus_b ( E + j * (d * d * k + p * d) + d * d * k, prod_ab (A + i * d * d + j * (d * d * k + d * p), d, d, E + i * (d * d * k + p * d) + d * d * k, d, p, TEMP), d, p);
		        }
        }
        
        pthread_barrier_wait(inform->barrier);//WAIT ALL
        
    }
    
																																																																																																																																											
    
    	double*IN=new double[p*p];
		double*T1=new double[p*p];
		 if(p!=0){
		 	double*TEMP=new double[p*d];
		 	copy(TEMP,A+k*(d*d*k+p*d)+p*d*k,p,p);
		    if (inverse_mat(TEMP, IN, p, *inform -> N)==0){
		        printf("ERROR in the last block\n");
		        pthread_mutex_lock(inform->error);
		        *inform->err=-3;
		        pthread_mutex_unlock(inform->error);
				delete[] TEMP;
				delete[] Ed;
				delete[] IN;
				delete[] T1;
				pthread_barrier_wait(inform->barrier);
		        return NULL;
		    }

		   //printf("INVERSE:\n");
		    //print_mat(IN,p);
		    //printf("\n");
		   //pthread_barrier_wait (inform -> barrier);
		   for (int j = k - 1 -id; j >=0 ; j-=t)
		   {                                         //prod E  last line to inverse element
		       prod_ab(IN,p,p,E+k*(d*d*k+p*d)+j*d*p,p,d,TEMP);
		       copy(E+k*(d*d*k+p*d)+j*p*d,TEMP,p,d);
		   }

		  	delete[] TEMP;
		  	if(id == 0)
		  	{
		    	prod_ab(IN,p,p,E+k*(d*d*k+p*d)+k*d*p,p,p,T1);
		    	copy(E+k*(d*d*k+p*d)+k*p*d,T1,p,p);
			}
		   	pthread_barrier_wait (inform -> barrier);
		   	
		    for(int i = k - 1 - id; i >= 0 ;i -= t){
		        for (int j=0;j<k;j++){
		            for(int i1=0;i1<d;i1++){
		                for(int i2=0;i2<d;i2++){
		                    for(int t1=0;t1<p;t1++){
		                        E[i*(d*d*k+d*p)+j*d*d+i1*d+i2]-=A[i*(d*d*k+d*p)+d*d*k+i1*p+t1]*E[k*(d*d*k+d*p)+j*d*p+d*t1+i2];
		                    }
		                }
		            }
		        }
		    }

		    for(int i = k - 1 -id; i >= 0; i -= t){

		        for(int i1=0;i1<d;i1++){
		            for(int i2=0;i2<p;i2++){
		                for(int t1=0;t1<p;t1++){
		                    E[i*(d*d*k+d*p)+k*d*d+i1*p+i2]-=A[i*(d*d*k+d*p)+d*d*k+i1*p+t1]*E[k*(d*d*k+d*p)+k*d*p+p*t1+i2];
		                }

		            }

		        }
		    }
		}
		
		
		//print_convert_mat(E,d,k,p);
		//printf("TUTUTUT\n");
		delete[] T1;
		delete[]IN;

    
    
    t0=get_time()-t0;
    t1=get_full_time()-t1;
    inform->TIME[id]=t0;
    pthread_barrier_wait(inform->barrier);
    if(fabs(*inform->err+3)<EPS){
				delete[] TEMP;
				delete[] Ed;
				//delete[] IN;
				//delete[] T1;
				//pthread_barrier_wait(inform->barrier);
		        return NULL;
    }
    /*delete[] TTT;
    delete[] TEMP;
    delete[] Ed;
    return NULL;*/
    if(id==0)
    {
    	printf( "ELAPSED: %f\n", t1);
    	for(int i = 0; i < t; i++)
    	{
    		printf( "T[%d] = %f\n",i ,inform->TIME[i] );
    	}
    	
    	CONVERTATION_STRING( inform->ind_col, k );
    	int n = k * d + p;
    	
    	if(k!=0)
    	{
        	convert_E(E,A,d,p,k,inform->ind_col);
        	copy(E,A,n,n);
    	}
    	if(argc==5){
		    f = fopen(argv[4], "r");
			printf("%s\n",argv[4]);
					
		    if(f == NULL){                                             //bad file
		        printf("CAN'T OPEN FILE\n");
				
		        
		        pthread_mutex_lock(inform->error);
		        	*inform->err=-2;
		        pthread_mutex_unlock(inform->error);
		        
		        delete[] TEMP;
		        delete[] Ed;
			    pthread_barrier_wait(inform->barrier);
		        return NULL;
		    }

		    if(f!=NULL && (read(A,f,k*d+p))<0){
		        printf("ERROR of reading or \n");
		        if(f!=NULL)
		        	fclose(f);
		        
		        pthread_mutex_lock(inform->error);
		        	*inform->err=-2;
		        pthread_mutex_unlock(inform->error);
		        delete[] TEMP;
		        delete[] Ed;
			    pthread_barrier_wait(inform->barrier);
			    return NULL;
		    }
		    if(f!=NULL)
		    	fclose(f);
		   // print_convert_mat(A,d,k,p);
		}
		else{
		    func_mat(A,k,d,p);
		}

		if( argc == 4)
		{
			 double *produ=new double[n*n];
			 if(produ==NULL){
				printf("CAN't convert\n");
				
				delete[] TEMP;
				delete[] Ed;			
				return 0;
			} 


			
				convert(A,produ,d,p,k);
				copy(A,produ,n,n);

			delete[] produ;
		}
		for(int i = 0; i < (k * d + p) * (k * d + p); i++){
			inform->res[i]=0.0;
		}
        if(n > 100 && t == 1){
            printf("without residual\n");
            pthread_mutex_lock(inform->error);
                *inform->err=-5;

            pthread_mutex_unlock(inform->error);
            delete[] TEMP;
            delete[] Ed;
            pthread_barrier_wait(inform->barrier);
            return NULL;

        }
		//print_mat(A,n);
		//print_mat(E,n);		
    }
    
    //printf("id = %d\n",id);
    pthread_barrier_wait(inform->barrier);
    if(( *inform->err + 2 ) < EPS){
        return NULL;
    }
    if(( *inform->err + 5 ) < EPS){
        return NULL;
    }
    
    
    inf_res inf;
        inf.A=A;
        inf.Inv=E;
        inf.t=t;
        inf.id=id;
        inf.res=inform->res;
        inf.n = d * k + p;
        inf.barrier=inform->barrier;
    
    count_residual(&inf);
    pthread_barrier_wait(inform->barrier);
    double *res1=inform->res;
    if(id==0){
    	int n=k*d+p;
		k=(d*k+p)/t;
		for(int i=t*k;i<n;i++){//peace we need to prod to the string
			for(int j1=0;j1<n;j1++){
				double s=0.0;
				for(int k1=0;k1<n;k1++){
					s+=A[i*n+k1]*E[k1*n+j1];
				}
				res1[i*n+j1]=s;
			}
		}
		
		for(int i=0;i<t*k;i++){
			for(int j1=t*k;j1<n;j1++){
				double s=0.0;
				for(int k1=0;k1<n;k1++){
					s+=A[i*n+k1]*E[k1*n+j1];
				}
				res1[i*n+j1]=s;
			}
			
		}
		//printf("\n");
		//print_mat(res,n);
		
		
		for(int i=0;i<n;i++){
		    res1[i*n+i]-=1;
		}
		//
		double r=norm(res1,n);
		
		printf("RESIDUAL: %e\n",r);
    
	}
    delete[] TEMP;
    delete[] Ed;
    
    return NULL;////
}




int search_block(double *A,double *E,int line,int k,int d,int p,int t,double N,int id){
    double Norm_min=0;
    bool ttt=0;
    double temp;
    int dex=0;
    double *T=new double[d*d];

    for (int i=k-1-id ; i>=line ; i-=t){//search not from begining
        copy(T,A+line*(d*d*k+d*p)+(d*d)*i,d,d);
        if(inverse_mat(T,E,d,N)){
            if(ttt==0){
                Norm_min=norm(E,d);
                dex=i;
                ttt=1;
                continue;

            }
            temp=norm(E,d);
            if(temp<Norm_min){
                Norm_min=temp;
                dex=i;
            }
        }
    }

    if(ttt==0){
        delete [] T;
        return -1;
    }
    copy(T,A+line*(d*d*k+p*d)+d*d*dex,d,d);
    inverse_mat(T,E,d,N);
    delete [] T;
    return dex;
}





int Jordan(double *A,double *E,int k,int p,int d,int *index_column,int t,pthread_t *id,int argc,char** argv){
    inform_thread* THR=new inform_thread[t];
    if(THR==NULL){
        printf("ERRROR memory\n");
        return -1;
    }
	
    double* aa=new double[d*d];
    if(aa==NULL){
        delete[] THR;
        printf("ERRROR memory\n");
        return -1;
    }
    pthread_mutex_t bb=PTHREAD_MUTEX_INITIALIZER;

    int* cc=new int[1];
    if(cc==NULL){
        delete[]THR;
        delete[]aa;
        printf("ERRROR memory\n");
        return -1;
    }
    double nno;
    int col;
    pthread_mutex_t dd=PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t ee;
    pthread_barrier_init(&ee,NULL,t);
    pthread_mutex_t ff=PTHREAD_MUTEX_INITIALIZER;
    nno=0;
    *cc=1;
    col=-1;
    double *N=new double[1];
    double * time=new double[t];
    int n=k*d+p;
    double *res= new double [n*n];
    if(res==NULL){
        delete[]THR;
        delete[] cc;
        delete[]aa;
        delete[] N;
        delete[] time;
        return -10;
    }
    for(int i=0;i<t;i++){
        THR[i].inverse_share=aa;
        THR[i].inverse_mutex=&bb;
        THR[i].err=cc;
        THR[i].error=&dd;
        THR[i].barrier=&ee;
        THR[i].id=i;
        THR[i].t=t;
        THR[i].res=res;
		THR[i].argc=argc;
		THR[i].TIME=time;
		THR[i].argv=argv;
		THR[i].N=N;
		
        THR[i].column=&col;
        THR[i].k=k;
        THR[i].p=p;
        THR[i].d=d;
        THR[i].norm=&nno;
        THR[i].help_mutex=&ff;
        THR[i].A=A;
        THR[i].E=E;

        THR[i].ind_col=index_column;
    }




    int temp;
    
	cpu_set_t cpu;
    //double t0=get_full_time();
    for(int i=0;i<t;i++){
        temp=pthread_create(id+i,NULL,&thread_0,THR+i);
        if(temp<0){
            printf("ERROR\n");
            delete[] N;
			delete[]time;
            pthread_barrier_destroy(&ee);
            delete[] THR;
            delete[] aa;
            delete[] cc;
            delete[] res;
            return -1;
        }
        CPU_ZERO(&cpu);
		CPU_SET(i,&cpu);
		pthread_setaffinity_np(id[i],sizeof(cpu),&cpu);
    }

    for(int i=0;i<t;i++){
        pthread_join(*(id+i),NULL);
        //printf("%d\n",i);
    }
    //t0=get_full_time()-t0;
	if(fabs(*THR[0].err+1.0)<EPS){
		printf("ERROR in the thread\n");

		delete[] THR;
		delete[] res;
		delete[]time;
        pthread_barrier_destroy(&ee);
		delete[] aa;
        delete[] cc;
        delete[] N;
		return -1;
	}
	
	if(fabs(*THR[0].err+2.0)<EPS){
        printf("ERROR in the file\n");
		delete[] N;
		delete[] res;
		delete[]time;
        pthread_barrier_destroy(&ee);
		delete[] THR;
		delete[] aa;
        delete[] cc;;
		printf("HAHAH\n");
		return -1;
	}
	if(fabs(*THR[0].err+3.0)<EPS){
		printf("ERROR in the last block\n");

		delete[] N;
		delete[] res;
		delete[]time;
        pthread_barrier_destroy(&ee);
		delete[] THR;
		delete[] aa;
        delete[] cc;
		return -1;
	}
	//printf("Elapsed: %f\n",t0);
	/*for(int i=0;i<t;i++){
		printf("ID=%d, Time: %f\n",i+1,time[i]);
	}*/




    




    //print_convert_mat(E,d,k,p);

    pthread_barrier_destroy(&ee);
    delete[] THR;
    delete[] aa;
    delete[] cc;
    delete[] N;
    delete[]time;
    delete[] res;
    return 1;
    //printf("%d",(int)id[0]);
}




int main(int argc, char** argv){

    if(argc<=3 || argc>=6){                                                //bad parametrs in command line
        printf("USE: ./a.out    size    block   threads  \nor\n./a.out    size    block   threads file\n");
        return -1;
    }


    int n;
    int d;
    int t;
    sscanf(argv[1],"%d",&n);
    sscanf(argv[2],"%d",&d);
    sscanf(argv[3],"%d",&t);



    double *mat=new double[n*n];
    for(int i=0;i<n*n;i++)
    	mat[i]=-10;
	if(mat==NULL){
        printf("ERROR of memory\n");
        return -4;
    }
    double *Ed=new double[n*n];
    if(Ed==NULL){
        printf("ERROR of memory\n");
        delete[]mat;
        return -4;
    }
     if(d <= 0 || n<=0 || t<=0){                                                  //bad block size
        printf("BAD parametrs\n");
        delete[] Ed;
        delete[]mat;
        return -3;
    }
    if(d%3!=0)
        d+=(3-d%3);
    if(n<=0){
        printf("BAD MATIRX SIZE\n");

    	delete[] Ed;
		delete[]mat;

        return -3;
    }
    int k=n/d,p=n%d;

    int *index_column=new int[k];
    if(index_column==NULL){
        printf("ERROR of memory\n");
        delete[] Ed;
        delete[]mat;
        return -2;
    }

    for(int i=0;i<k;i++)
        index_column[i]=i;
	
    //printf("%d\n",d);



    
    pthread_t *id=new pthread_t[t];
    if(id==NULL){
        printf("Error memory\n");
        delete[] mat;
        delete[] Ed;
        delete[] index_column;
        return -3;
    }

    if(Jordan(mat,Ed,k,p,d,index_column,t,id,argc,argv)<0){
        printf("JORDAN CAN't WORK\n");
        delete[] mat;
        delete[] id;
        delete[] Ed;
        delete[] index_column;
        return -3;
    }
    
    
  // PRINT_INVERSE(Ed,d,k,p,index_column);
	
/*    
    printf("INVERSE MATRIX:\n");
    print_mat(Ed,n);
    printf("\nINITIAL MATRIX:\n");
    print_mat(mat,n);
	*/
    /*prod_ab(mat,n,n,Ed,n,n,produ);
    for(int i=0;i<n;i++){
        produ[i*n+i]-=1.0;
    }
    printf("\nRESULT:\n");
    //print_mat(produ,n);
    double tt=norm(produ,n);
    prod_ab(Ed,n,n,mat,n,n,produ);
    double tt1=norm(produ,n);
    if(tt>tt1)
        tt=tt1;
    printf("RESIDUAL IN MAIN: %f\n",tt);
    */
    //residual(mat,Ed,n,t);
    //sleep(2);
    
    delete[] id;
    delete[] Ed;
    delete []index_column;
    delete []mat;
   //fclose(f);

    return 0;
}
































































































int typical_prod(double *a,int a1,int a2, double *b,int b1,int b2,double *c){
    for(int i=0;i<a1;i++){
        for(int j=0;j<b2;j++){
            double s=0.0;
            for(int k=0;k<a2;k++){
                s+=a[i*a2+k] * (b[k*b2+j]);
            }
            c[i*b2+j]=s;
        }
    }
    return 0;
    printf("%d",b1);
}





void CONVERTATION_STRING(int *ind,int n){
    int *t=new int[n];
    for(int i=0;i<n;i++){
        t[ind[i]]=i;
    }
    for(int i=0;i<n;i++){
        ind[i]=t[i];
    }
    delete[] t;

}

int read(double* mat, FILE *f,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(fscanf(f,"%lf",&mat[i*n+j])!=1){
                printf("ERROR of fscanf\n");
                return -1;
            }
        }
    }
    return 1;
}

void PRINT_INVERSE(double *a,int d,int k,int p,int *ind){
    for(int j=0;j<k;j++){
        for(int i=0;i<k;i++){

            for(int j1=0;j1<d;j1++){
                for(int j2=0;j2<d;j2++)
                    printf("%f ",a[ind[j]*(d*d*k+p*d)+d*d*i+j1*d+j2]);
                printf("\n");
            }
            printf("\n");
        }
         for(int j1=0;j1<d;j1++){
                for(int j2=0;j2<p;j2++)
                    printf("%f ",a[ind[j]*(d*d*k+p*d)+d*d*k+j1*p+j2]);
            printf("\n");

        }
        printf("\n\n");
    }
    for(int j=0;j<k;j++){
        for(int j1=0;j1<p;j1++){
            for(int j2=0;j2<d;j2++)
                printf("%f ",a[k*(d*d*k+p*d)+d*p*j+j1*d+j2]);
            printf("\n");
        }
        printf("\n");
    }
    for(int j1=0;j1<p;j1++){
        for(int j2=0;j2<p;j2++){
            printf("%f ",a[k*(d*d*k+p*d)+d*p*k+j1*p+j2]);
        }
        printf("\n");
    }
    printf("***********************\n***********************\n");

}




void func_mat(double *A,int k, int d, int p){
    for(int i=0;i<k;i++){
        for(int j=0;j<k;j++){
            for(int i1=0;i1<d;i1++){
                for(int i2=0;i2<d;i2++){
                    A[i*(d*d*k+d*p)+j*d*d+i1*d+i2]=fabs(i*d+i1-(j*d+i2));//1.0/(1+(i*d+i1+j*d+i2));
                }
            }
        }
    }
    for(int i=0;i<k;i++){
        for(int i1=0;i1<d;i1++){
            for(int i2=0;i2<p;i2++){
                    A[i*(d*d*k+d*p)+k*d*d+i1*p+i2]=fabs(i*d+i1-(k*d+i2));//1.0/(1+(i*d+i1+k*d+i2));
            }
        }
    }
    for(int i=0;i<k;i++){
        for(int i1=0;i1<p;i1++){
            for(int i2=0;i2<d;i2++){
                    A[k*(d*d*k+d*p)+i*d*p+i1*d+i2]=fabs(k*d+i1-(i*d+i2));//1.0/(1+(k*d+i1+i*d+i2));
            }
        }
    }
    for(int i1=0;i1<p;i1++){
        for(int i2=0;i2<p;i2++){
            A[k*(d*d*k+d*p)+k*d*p+i1*p+i2]=fabs(k*d+i1-(k*d+i2));//1.0/(1+(k*d+i1+k*d+i2));
        }
    }
}
/*
void func_mat(double *A,int k, int d, int p){
    for(int i=0;i<k;i++){
        for(int j=0;j<k;j++){
            for(int i1=0;i1<d;i1++){
                for(int i2=0;i2<d;i2++){
                    A[i*(d*d*k+d*p)+j*d*d+i1*d+i2]=1.0/(1+(i*d+i1+j*d+i2));
                }
            }
        }
    }
    for(int i=0;i<k;i++){
        for(int i1=0;i1<d;i1++){
            for(int i2=0;i2<p;i2++){
                    A[i*(d*d*k+d*p)+k*d*d+i1*p+i2]=1.0/(1+(i*d+i1+k*d+i2));
            }
        }
    }
    for(int i=0;i<k;i++){
        for(int i1=0;i1<p;i1++){
            for(int i2=0;i2<d;i2++){
                    A[k*(d*d*k+d*p)+i*d*p+i1*d+i2]=1.0/(1+(k*d+i1+i*d+i2));
            }
        }
    }
    for(int i1=0;i1<p;i1++){
        for(int i2=0;i2<p;i2++){
            A[k*(d*d*k+d*p)+k*d*p+i1*p+i2]=1.0/(1+(k*d+i1+k*d+i2));
        }
    }
}
*/




void convert_E(double *a,double *b,int d,int p,int k,int *ind){
    for(int i=0;i<k;i++){//middle
        for(int j=0;j<k;j++){
            for(int i1=0;i1<d;i1++){
                for(int i2=0;i2<d;i2++){
                    b[(i*d+i1)*(k*d+p)+d*j+i2]=a[ind[i]*(d*d*k+p*d)+j*d*d+i1*d+i2];
                }
            }
        }
    }
    for(int i=0;i<k;i++){//right
         for(int i1=0;i1<d;i1++){
                for(int i2=0;i2<p;i2++){
                        b[(i*d+i1)*(k*d+p)+d*k+i2]=a[ind[i]*(d*d*k+p*d)+k*d*d+i1*p+i2];
                }

        }
    }
    for(int i=0;i<k;i++){//down
         for(int i1=0;i1<p;i1++){
                for(int i2=0;i2<d;i2++){
                        b[(k*d+i1)*(k*d+p)+d*i+i2]=a[k*(d*d*k+p*d)+i*d*p+i1*d+i2];
                }

        }
    }

    for(int i1=0;i1<p;i1++){
                for(int i2=0;i2<p;i2++){
                        b[(k*d+i1)*(k*d+p)+d*k+i2]=a[k*(d*d*k+p*d)+k*d*p+i1*p+i2];
                }

        }

}

void convert(double *a,double *b,int d,int p,int k){
    for(int i=0;i<k;i++){//middle
        for(int j=0;j<k;j++){
            for(int i1=0;i1<d;i1++){
                for(int i2=0;i2<d;i2++){
                    b[(i*d+i1)*(k*d+p)+d*j+i2]=a[i*(d*d*k+p*d)+j*d*d+i1*d+i2];
                }
            }
        }
    }
    for(int i=0;i<k;i++){//right
         for(int i1=0;i1<d;i1++){
                for(int i2=0;i2<p;i2++){
                        b[(i*d+i1)*(k*d+p)+d*k+i2]=a[i*(d*d*k+p*d)+k*d*d+i1*p+i2];
                }

        }
    }
    for(int i=0;i<k;i++){//down
         for(int i1=0;i1<p;i1++){
                for(int i2=0;i2<d;i2++){
                        b[(k*d+i1)*(k*d+p)+d*i+i2]=a[k*(d*d*k+p*d)+i*d*p+i1*d+i2];
                }

        }
    }

    for(int i1=0;i1<p;i1++){
                for(int i2=0;i2<p;i2++){
                        b[(k*d+i1)*(k*d+p)+d*k+i2]=a[k*(d*d*k+p*d)+k*d*p+i1*p+i2];
                }

        }

}

void print_convert_mat(double *a,int d,int k,int p){
    for(int j=0;j<k;j++){
        for(int i=0;i<k;i++){

            for(int j1=0;j1<d;j1++){
                for(int j2=0;j2<d;j2++)
                    printf("%10.10f ",a[j*(d*d*k+p*d)+d*d*i+j1*d+j2]);
                printf("\n");
            }
            printf("\n");
        }
         for(int j1=0;j1<d;j1++){
                for(int j2=0;j2<p;j2++)
                    printf("%10.10f ",a[j*(d*d*k+p*d)+d*d*k+j1*p+j2]);
            printf("\n");

        }
        printf("\n\n");
    }
    for(int j=0;j<k;j++){
        for(int j1=0;j1<p;j1++){
            for(int j2=0;j2<d;j2++)
                printf("%10.10f ",a[k*(d*d*k+p*d)+d*p*j+j1*d+j2]);
            printf("\n");
        }
        printf("\n");
    }
    for(int j1=0;j1<p;j1++){
        for(int j2=0;j2<p;j2++){
            printf("%10.10f ",a[k*(d*d*k+p*d)+d*p*k+j1*p+j2]);
        }
        printf("\n");
    }
    printf("***********************\n***********************\n");

}


void null(double *a,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            a[i*n+j]=0.0;
        }
    }

}

void unit(double *a,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            a[i*n+j]=(double)(i==j);
        }
    }

}

void copy(double *a,double *b,int n,int m){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)
            a[i*m+j]=b[i*m+j];
    }
}




void a_minus_b(double *a,double *b,int n,int m){
    for(int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            a[i*m+j]-=b[i*m+j];
        }
    }
}

int inverse_mat(double *A, double* E, int m,double N){
    int i,j,k,max_i;
    double buf, max,tmp;


    for(i=0;i<m;i++)
        for(j=0;j<m;j++)
            E[i*m+j]=(double)(i==j);

    for(i=0;i<m;i++){
        max=fabs(A[i*m+i]);
        max_i=i;
        for(j=i+1;j<m;j++){//max in column+index
            if(fabs(A[j*m+i])>max){
                            max=fabs(A[j*m+i]);
                            max_i=j;
            }
        }
        if(max<EPS*N)
            return 0;

        if(max_i!=i){//swap a part of line in both matrix
            for(k=i;k<m;k++){
                buf=A[i*m+k];
                A[i*m+k]=A[max_i*m+k];
                A[max_i*m+k]=buf;
                buf=E[i*m+k];
                E[i*m+k]=E[max_i*m+k];
                E[max_i*m+k]=buf;
            }
            for(k=0;k<i;k++){//swap another part of line in E matrix
                buf=E[i*m+k];
                E[i*m+k]=E[max_i*m+k];
                E[max_i*m+k]=buf;
            }
        }

        buf=1.0/A[i*m+i];
        for(j=i+1;j<m;j++)
            A[i*m+j]*=buf;

        for(j=0;j<m;j++)
            E[i*m+j]*=buf;

        for(j=i+1;j<m;j++){
            buf=A[j*m+i];
            for(k=i+1;k<m;k++)
                A[j*m+k]-=buf*A[i*m+k];
            for(k=0;k<m;k++)
                E[j*m+k]-=buf*E[i*m+k];
        }
    }

    for (k=0;k<m;k++){
        for (i=m-1;i>=0;i--){
            tmp=E[i*m+k];
            for(j=i+1;j<m;++j)
                tmp-=A[i*m+j]*E[j*m+k];
            E[i*m+k]=tmp;
        }
    }
    return 1;
}



double* prod_ab(double* a,int a1,int a2, double *b,int b1,int b2,double*c){

     double s[9];

    for(int i=1;i<a1/3 *3;i+=3){
        for(int j=1;j<b2/3 *3;j+=3){
            for(int k=0;k<9;k++)
                s[k]=0;
            for(int k=0;k<a2;k++){
                s[0]+=a[(i-1)*a2 + k]*b[(j-1)+b2*k];
                //printf("%f %f \n",a[(i-1)*a2 + k],b[(j-1)+b2*k]);
                s[1]+=a[(i-1)*a2 + k]*b[(j)+b2*k];
                s[2]+=a[(i-1)*a2 + k]*b[(j+1)+b2*k];
                s[3]+=a[(i)*a2 + k]*b[(j-1)+b2*k];
                s[4]+=a[(i)*a2 + k]*b[(j)+b2*k];
                s[5]+=a[(i)*a2 + k]*b[(j+1)+b2*k];
                s[6]+=a[(i+1)*a2 + k]*b[(j-1)+b2*k];
                s[7]+=a[(i+1)*a2 + k]*b[(j)+b2*k];
                s[8]+=a[(i+1)*a2 + k]*b[(j+1)+b2*k];

            }
           // printf("\n");
            c[(i-1)*b2+(j-1)]=s[0];
            c[(i-1)*b2+(j)]=s[1];
            c[(i-1)*b2+(j+1)]=s[2];
            c[(i)*b2+(j-1)]=s[3];
            c[(i)*b2+(j)]=s[4];
            c[(i)*b2+(j+1)]=s[5];
            c[(i+1)*b2+(j-1)]=s[6];
            c[(i+1)*b2+(j)]=s[7];
            c[(i+1)*b2+(j+1)]=s[8];
        }
    }

    for(int i=0;i<(a1/3)*3;i++){
        for(int j=(b2/3)*3;j<b2;j++){
            s[0]=0;
            for(int k=0;k<a2;k++){
                s[0]+=a[i*a2+k]*b[k*b2+j];
            }
            c[i*b2+j]=s[0];
        }
    }

    for(int i=(a1/3)*3;i<a1;i++){
        for(int j=0;j<b2;j++){
            s[0]=0;
            for(int k=0;k<a2;k++){
                //printf("%f %f \n",a[i*a2+k],b[k*b2+j]);
                s[0]+=a[i*a2+k]*b[k*b2+j];
            }
            c[i*b2+j]=s[0];

        }
        //printf("\n");

    }

    return c;
    printf("%d",b1);
}


double norm(double *a,int n){
    double s=0,t=0;
    for(int i=0;i<n;i++){
        s=0;
        for(int j=0;j<n;j++){
            s+=fabs(a[i*n+j]);
        }
        if(s>t)
            t=s;
    }
    return t;

}


int read(double* mat,FILE *f,int n,int d){

   // printf("%d azaza \n",*n);
    double temp;
    int n_bl=(n) / d;
    int p=n%d;
   // printf("%d \n\n",n_bl);
    for(int i = 1; i <= (n); i++){
        for(int j = 1; j <= (n); j++){
            if(fscanf(f,"%lf",&temp)!=1){
                printf("ERROR, CAN'T READ MATRIX ELEMENTS\n");
            
                return -1;
            }
			
        	//printf("i=%d j=%d\n%f\n",i,j,temp);
            int n1,n2;
            int k1=i/d,k2=i%d;
            if(!(k1 ==n_bl && ( k2 != 0))){
              n1=k1;
              if(k2 != 0)
                  n1++;
            }
            else
                n1=n_bl+1;

            k1=j/d;
            k2=j%d;
            if(!( k1 ==n_bl && ( k2 != 0))){
              n2=(j)/d;
              if(k2 != 0)
                  n2++;
            }
            else
                n2=n_bl+1;

            //printf("%d %d \n",n1,n2);

            if(n1<=n_bl){
                if(n2>n_bl){//right block
                    mat[(n1-1)*(d*d*n_bl+d*p)+(n2-1)*(d*d)+((i-1)%d)*p+((j-1)%d)]=temp;

                }
                else{//normal block
                    mat[(n1-1)*(d*d*n_bl+d*p)+(n2-1)*(d*d)+((i-1)%d)*d+((j-1)%d)]=temp;
                }
            }
            else{
                if(n2<=n_bl){//low block
                    mat[(n1-1)*(d*d*n_bl+d*p)+(n2-1)*(p*d)+((i-1)%d)*d+((j-1)%d)  ]=temp;
                }
                else{//low right block
                     mat[(n1-1)*(d*d*n_bl+d*p)+(n2-1)*(p*d)+((i-1)%d)*p+((j-1)%d)  ]=temp;
                }
            }

        }

    }
    return 1;

}
void print_mat(double *a,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            printf("%10.10f ",a[i*n+j]);
        }
        printf("\n");

    }

}

