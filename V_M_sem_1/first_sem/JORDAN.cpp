#include<pthread.h>
#include<stdio.h>
#include<time.h>
#include<cmath>
#include<sys/time.h>
#include <unistd.h>

#define EPS 1e-10
//in.txt n d     or      n d
//-ffast-math -O3 -march=native
double* read(FILE *f,int n,int d);//OK
void print_mat(double *a,int n);//OK
void print_convert_mat(double *a,int d,int k,int p);//OK
double* prod_ab(double* a,int a1,int a2, double *b,int b1,int b2,double* c);//OK
double norm(double *a,int n);//OK
int inverse_mat(double* A,double* E,int m,double N);//OK inverse matrix in E      //1-Ok, 0-Err
void a_minus_b(double *a,double *b,int n,int m);//OK
int Jordan(double *A,double *E,int k,int p,int d,int *index_column,double N, int t, pthread_t *id);

int search_block(double *A,double *E,int line,int k,int d,int p,int t,double N,int id);
void copy(double *a,double *b,int n,int m);//OK
void null(double *a,int n);//OK
void unit(double *a,int n);//OK
void convert(double *a,double *b,int d,int p,int k);//convertation from  a in b
void func_mat(double *A,int k, int d, int p);//GILBERT
void PRINT_INVERSE(double *a,int d,int k,int p,int *ind);//PRINT AFTER TRANSFORMATION
double* read(FILE *f,int n);
void convert_E(double *a,double *b,int d,int p,int k,int *ind);
void CONVERTATION_STRING(int *ind,int n);
void* thread_0(void* T);
int typical_prod(double *a,int a1,int a2, double *b,int b1,int b2,double *c);
void Jord(void*T);
struct inform_thread{
    double *A;
    double *E;
    int* ind_col;
    double N;
    int p;
    int k;
    int t;
    int d;

    int *column;
    double* inverse_share;
    double* norm;
    pthread_mutex_t *inverse_mutex;

    int* err;
    pthread_mutex_t *error;

    pthread_barrier_t *barrier;
    int id;
    pthread_mutex_t* help_mutex;
};


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

void Jord(void*T){
	inform_thread *inform=(inform_thread*)(T);
    int d=inform->d;
    int p=inform->p;
    int k=inform->k;
    //int n=k*d+p;
    int t=inform->t;
    int id=inform->id;
    double* Ed=new double[d*d];
    double* TEMP=new double[d*d];
    double* E=inform->E;
    double *A=inform->A;
    int res=0;
    printf("id=%d \n",id);
    double *TTT=new double[d*d*k+p*d];
   // double time=get_cpu_time();
    for(int i=0;i<k;i++){

        res=search_block(A,Ed,i, k, d, p, t, inform->N, id);//returns -1 if cannot find inverse
       /*	
       	////////////////////////////////////////////////////////////////////////PPPPPPPPPPRRRRRRRRRRRIIIIIIIIIIIINNNNNNNNNNT
       	pthread_mutex_lock(inform->inverse_mutex);
        printf("ЕГЕГЕГЕГЕГ\n***********************************************************************\nid=%d\n",id);
        printf("INVERSE:\n");
        printf("\n");
        print_mat(Ed,d);
        pthread_mutex_unlock(inform->inverse_mutex);
        
       	////////////////////////////////////////////////////////////////////////PPPPPPPPPPRRRRRRRRRRRIIIIIIIIIIIINNNNNNNNNNT
       	*/ 
       	
       	
       if(res>=0){
        	
        
            pthread_mutex_lock(inform->inverse_mutex);
            if(fabs(*(inform->norm))<EPS){
                *(inform->column)=res;
                //printf("HAHAH\n\n\n");
                *(inform->norm)=norm(Ed,d);
                //printf("HAHAH %lf\n\n\n",*(inform->norm));
            }
            else{
                double nn=norm(Ed,d);
                if(*(inform->norm)>nn){
                    *(inform->column)=res;
                   // printf("HAHAH %lf\n\n\n",nn);
                    *(inform->norm)=nn;
                }
            }


            pthread_mutex_unlock(inform->inverse_mutex);

        }
        
        pthread_barrier_wait(inform->barrier);//WAIT INVERSE
        
        if( id == 0){
		    
		    if(!(fabs(*inform->norm) <  EPS) ){//minus
		        //printf("SWAPSWAPSWAP\n\n");
		        *inform->norm=0.0;//i and res
		       // printf("i=%d res=%d\n",i,res);
		        
		        //printf("ind[i]=%d ind[res]=%d\n",inform->ind_col[i],inform->ind_col[res]);
		        
		        int ttt=inform->ind_col[i];
		        inform->ind_col[i]=inform->ind_col[*inform->column];
		        inform->ind_col[*inform->column]=ttt;

		    }
		    else{
		        *inform->err=-1;
		        *(inform->column)=1;
		    }
        }
        
		
       
        for(int j=k-1-id;j>=0;j-=t){
            copy(TEMP,A+j*(d*d*k+p*d)+*inform->column*d*d,d,d);
            copy(A+j*(d*d*k+p*d)+*inform->column*d*d,A+j*(d*d*k+p*d)+i*d*d,d,d);
            copy(A+j*(d*d*k+p*d)+i*d*d,TEMP,d,d);
        }
        
        if(p!=0 && inform->id==0){
		    copy(TEMP,A+k*(d*d*k+p*d)+*inform->column*d*p,p,d);
		    copy(A+k*(d*d*k+p*d)+*inform->column*p*d,A+k*(d*d*k+p*d)+i*p*d,p,d);
		    copy(A+k*(d*d*k+p*d)+i*p*d,TEMP,p,d);
		}
		
        //pthread_barrier_wait(inform->barrier);//WAIT for copy
        
       	if(*inform->err==-1){
            delete[] TEMP;
            delete[] Ed;
            //return NULL;
        }
		//pthread_mutex_lock(inform->help_mutex);
        copy(TEMP,A+i*(d*d*k+p*d)+i*d*d,d,d);
        inverse_mat(TEMP,Ed,d,inform->N);
       
		//pthread_mutex_unlock(inform->help_mutex);
        
        /*
		////////////////////////////////////////////////////////////////////////PPPPPPPPPPRRRRRRRRRRRIIIIIIIIIIIINNNNNNNNNNT
       	if(id==0){
		   	pthread_mutex_lock(inform->inverse_mutex);
		    printf("FINAL INVERSE:\n");
		    printf("\n");
		    print_mat(Ed,d);
		    printf("\n");
		    printf("AFTER swap:\n");
		    print_convert_mat(A,d,k,p);
		    pthread_mutex_unlock(inform->inverse_mutex);
        }
       	////////////////////////////////////////////////////////////////////////PPPPPPPPPPRRRRRRRRRRRIIIIIIIIIIIINNNNNNNNNNT
        */ 
        for(int j=k-1-id ; j>i ; j-=t){                                         //prod line to inverse element
            prod_ab(Ed,d,d,A+i*(d*d*k+p*d)+j*d*d,d,d,TEMP);
            copy(A+i*(d*d*k+p*d)+j*d*d,TEMP,d,d);
        }


        for(int j=k-1-id ; j>=0 ; j-=t){                                         //prod E line to inverse element
            prod_ab(Ed,d,d,E+i*(d*d*k+p*d)+j*d*d,d,d,TEMP);
            copy(E+i*(d*d*k+p*d)+j*d*d,TEMP,d,d);
        }

        if(p!=0 && id==0){
            prod_ab(Ed,d,d,A+i*(d*d*k+p*d)+k*d*d,d,p,TEMP);
            copy(A+i*(d*d*k+p*d)+k*d*d,TEMP,d,p);

            prod_ab(Ed,d,d,E+i*(d*d*k+p*d)+k*d*d,d,p,TEMP);
            copy(E+i*(d*d*k+p*d)+k*d*d,TEMP,d,p);

        }
        
        
        
        pthread_mutex_lock(inform->help_mutex);
        for(int j=0;j<k;j++){
        	copy(TTT+j*d*d,A+i*d*d+j*(d*d*k+d*p),d,d);
        }
        	copy(TTT+k*d*d,A+i*d*p+k*(d*d*k+d*p),d,p);
        
        pthread_mutex_unlock(inform->help_mutex);
        
        
        for(int j=0;j<i;j++){// high correction
            for(int y=k-1-id;y>i;y-=t){

                a_minus_b(A+j*(d*d*k+p*d)+d*d*y,prod_ab(TTT+j*d*d,d,d,A+i*(d*d*k+p*d)+d*d*y,d,d,TEMP),d,d);
            }
        }////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for(int j=0;j<i;j++){// high E correction
            for(int y=k-1-id;y>=0;y-=t){

                a_minus_b(E+j*(d*d*k+p*d)+d*d*y,prod_ab(TTT+j*d*d,d,d,E+i*(d*d*k+p*d)+d*d*y,d,d,TEMP),d,d);
            }
        }

        for(int j=i+1;j<k;j++){//low correction
            for(int y=k-1-id;y>i;y-=t){

                a_minus_b(A+j*(d*d*k+p*d)+d*d*y,prod_ab(TTT+j*d*d,d,d,A+i*(d*d*k+p*d)+d*d*y,d,d,TEMP),d,d);
            }
        }

        for(int j=i+1;j<k;j++){//low E correction
            for(int y=k-1-id;y>=0;y-=t){

                a_minus_b(E+j*(d*d*k+p*d)+d*d*y,prod_ab(TTT+j*d*d,d,d,E+i*(d*d*k+p*d)+d*d*y,d,d,TEMP),d,d);
            }
        }

		
        if( p!=0){
            //printf("HAHAHA\n");
            for(int j=k-1-id;j>i;j-=t){//lowest correction
                    a_minus_b(A+k*(d*d*k+p*d)+d*p*j,prod_ab(TTT+k*d*d,p,d,A+i*(d*d*k+p*d)+d*d*j,d,d,TEMP),p,d);
            }
            if(id==0){
				for(int j=i+1;j<k;j++){//right correction
		                a_minus_b(A+j*(d*d*k+p*d)+d*d*k,prod_ab(TTT+j*d*d,d,d,A+i*(d*d*k+p*d)+d*d*k,d,p,TEMP),d,p);
		        }
            
			
		        for(int j=0;j<i;j++){//right correction
		                a_minus_b(A+j*(d*d*k+p*d)+d*d*k,prod_ab(TTT+j*d*d,d,d,A+i*(d*d*k+p*d)+d*d*k,d,p,TEMP),d,p);
		        }
		    
            	a_minus_b(A+k*(d*d*k+p*d)+d*p*k,prod_ab(TTT+k*d*d,p,d,A+i*(d*d*k+p*d)+d*d*k,d,p,TEMP),p,p);
            }
            for(int j=k-1-id;j>=0;j-=t){//lowest E correction
                    a_minus_b(E+k*(d*d*k+p*d)+d*p*j,prod_ab(TTT+k*d*d,p,d,E+i*(d*d*k+p*d)+d*d*j,d,d,TEMP),p,d);
            }
            if(id==0){
				for(int j=i+1;j<k;j++){//right E correction
						a_minus_b(E+j*(d*d*k+p*d)+d*d*k,prod_ab(TTT+j*d*d,d,d,E+i*(d*d*k+p*d)+d*d*k,d,p,TEMP),d,p);
				}
			   
				for(int j=0;j<i;j++){//right E correction
						a_minus_b(E+j*(d*d*k+p*d)+d*d*k,prod_ab(TTT+j*d*d,d,d,E+i*(d*d*k+p*d)+d*d*k,d,p,TEMP),d,p);
				}
		
			   	a_minus_b(E+k*(d*d*k+p*d)+d*p*k,prod_ab(TTT+k*d*d,p,d,E+i*(d*d*k+p*d)+d*d*k,d,p,TEMP),p,p);
			}
			

        }
        
        pthread_barrier_wait(inform->barrier);//WAIT ALL
        
    }
    
    delete[] TEMP;
    delete[] Ed;
}
void* thread_0(void* T){
    inform_thread* aa=(inform_thread*)T;
    struct timeval tm;
	double t0;
	gettimeofday(&tm,0);
	t0 = (double)tm.tv_sec + (double)tm.tv_usec * 1e-6 ;
    Jord(T);
    t0 = (double)tm.tv_sec + (double)tm.tv_usec * 1e-6 - t0;
	if (aa->id == 0) printf (" Elapsed = %.5f\n", t0); 
    return NULL;////
}









int Jordan(double *A,double *E,int k,int p,int d,int *index_column,double N,int t,pthread_t *id){
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
    int *col=new int[1];
    if(col==NULL){
        delete[]THR;
        delete[] cc;
        delete[]aa;
        printf("ERRROR memory\n");
        return -1;
    }
    double* nno=new double[1];
    if(nno==NULL){
        delete[]col;
        delete[]THR;
        delete[] cc;
        delete[]aa;
        printf("ERRROR memory\n");
        return -1;
    }
    pthread_mutex_t dd=PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t ee;
    pthread_barrier_init(&ee,NULL,t);
    pthread_mutex_t ff=PTHREAD_MUTEX_INITIALIZER;
    *nno=0;
    *cc=1;
    *col=-1;
    for(int i=0;i<t;i++){
        THR[i].inverse_share=aa;
        THR[i].inverse_mutex=&bb;
        THR[i].err=cc;
        THR[i].error=&dd;
        THR[i].barrier=&ee;
        THR[i].id=i;
        THR[i].t=t;

        THR[i].column=col;
        THR[i].k=k;
        THR[i].p=p;
        THR[i].d=d;
        THR[i].norm=nno;
        THR[i].help_mutex=&ff;
        THR[i].A=A;
        THR[i].E=E;

        THR[i].N=N;

        THR[i].ind_col=index_column;
    }




    int temp;
    

    for(int i=0;i<t;i++){
        temp=pthread_create(id+i,NULL,&thread_0,THR+i);
        if(temp<0){
            printf("ERROR\n");

            delete[]nno;
            delete[] THR;
            delete[] aa;
            delete[] cc;
            delete[]nno;
            delete[]col;
            return -1;
        }
    }

    for(int i=0;i<t;i++){
        pthread_join(*(id+i),NULL);
        //printf("%d\n",i);
    }
	if(fabs(*THR[0].err+1.0)<EPS){
		printf("ERROR in the thread\n");
		delete[]nno;
		delete[] THR;
		delete[] aa;
		delete[] cc;
		delete[]nno;
		delete[]col;
		return -1;
	}





    double*IN=new double[p*p];
    double*T1=new double[p*p];
     if(p!=0){
        if (inverse_mat(A+k*(d*d*k+p*d)+p*d*k,IN,p,N)==0){
            printf("ERROR in the last block\n");
            delete []IN;
            delete[] T1;
            return -2;
        }

       //printf("INVERSE:\n");
        //print_mat(IN,p);
        //printf("\n");
       double*TEMP=new double[p*d];
       for(int j=0;j<k;j++){                                         //prod E  last line to inverse element
           typical_prod(IN,p,p,E+k*(d*d*k+p*d)+j*d*p,p,d,TEMP);
           copy(E+k*(d*d*k+p*d)+j*p*d,TEMP,p,d);
       }

       delete[] TEMP;
        typical_prod(IN,p,p,E+k*(d*d*k+p*d)+k*d*p,p,p,T1);
        copy(E+k*(d*d*k+p*d)+k*p*d,T1,p,p);

       //printf("INVERSE:\n");
       //print_convert_mat(E,d,k,p);
     //   printf("\nHERE\n");
    //print_convert_mat(E,d,k,p,index_column);


        for(int i=0;i<k;i++){
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

        for(int i=0;i<k;i++){

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





    //print_convert_mat(E,d,k,p);

    pthread_barrier_destroy(&ee);
    delete[] THR;
    delete[] aa;
    delete[] cc;
    delete[] nno;
    delete[]col;
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


    FILE *f;
    double *mat;

    double *Ed=new double[n*n];
    if(Ed==NULL){
        printf("ERROR of memory\n");
        return -4;
    }
     if(d <= 0 || n<=0 || t<=0){                                                  //bad block size
        printf("BAD parametrs\n");
        delete[] Ed;
        return -3;
    }
    if(d%3!=0)
        d+=(3-d%3);
    if(n<=0){
        printf("BAD MATIRX SIZE\n");

    delete[] Ed;


        return -3;
    }
    int k=n/d,p=n%d;

    int *index_column=new int[k];
    if(index_column==NULL){
        printf("ERROR of memory\n");
        delete[] Ed;
        return -2;
    }

    for(int i=0;i<k;i++)
        index_column[i]=i;

    if(argc==5){
        f = fopen(argv[4], "r");

        if(f == NULL){                                             //bad file
            printf("CAN'T OPEN FILE\n");

        delete[] Ed;
        delete []index_column;

            return -2;
        }

        if((mat=read(f,n,d))==NULL){
            printf("ERROR of reading or memory\n");
            fclose(f);
        delete[] Ed;
        delete []index_column;

            return -3;
        }
        fclose(f);
    }
    else{
        mat=new double[n*n];
        func_mat(mat,k,d,p);
    }
    //printf("%d\n",d);



     double N=norm(mat,n);
    null(Ed,d*k+p);
    for(int i=0;i<k;i++){
        for(int j1=0;j1<d;j1++){

            Ed[i*(d*d*k+p*d)+i*d*d+j1*d+j1]=1.0;

        }
    }
    for(int i=0;i<p;i++){
        Ed[k*(d*d*k+p*d)+d*p*k+i*p+i]=1.0;
    }
    pthread_t *id=new pthread_t[t];
    if(id==NULL){
        printf("Error memory\n");
        delete[] mat;
        delete[] Ed;
        delete[] index_column;
        return -3;
    }
    unsigned int beg=clock();

    if(Jordan(mat,Ed,k,p,d,index_column,N,t,id)<0){
        printf("JORDAN CAN't WORK\n");
        delete[] mat;
        delete[] id;
        delete[] Ed;
        delete[] index_column;
        return -3;
    }
    //for(int i=0;i<k;i++)
    //    printf("%d",index_column[i]);
    unsigned int end=clock();
    printf("\n\n TIME: %f \n",(double)((end-beg))/CLOCKS_PER_SEC);
   // sleep(10);
    //printf("\n");
    //for(int i=0;i<k;i++){
     //   printf("%d ",index_column[i]);
   // }
   // printf("\n");
    CONVERTATION_STRING(index_column,k);
  // PRINT_INVERSE(Ed,d,k,p,index_column);

    if(k!=0){
        convert_E(Ed,mat,d,p,k,index_column);
        copy(Ed,mat,n,n);
   }
    beg=end;
    end=clock();
    //printf("\n\n TIME: %f \n",(double)((end-beg))/CLOCKS_PER_SEC);
    if(argc==5){


        f = fopen(argv[4], "r");
        if(f==NULL){
            printf("Can't open file twise\n");
            delete[] id;
            delete[] Ed;
            delete []index_column;
            delete []mat;

            return -2;
        }

        delete []mat;
        mat=read(f,n);
        if(f==NULL)
            printf("HUIHUIHUI\n");
        if(mat==NULL){
            printf("Can't read twise\n");
            delete[] id;
            delete[] Ed;
            delete []index_column;
            fclose(f);
            return -2;

        }
         fclose(f);

     //printf("HAHAHAH\n");
    }
    else{
        func_mat(mat,k,d,p);
    }
     //printf("HAHAHAH\n");
    double *produ=new double[n*n];
     if(produ==NULL){
        delete[] id;
        delete[] Ed;
        delete []index_column;
        delete []mat;
        return 0;
    }


    if(argc==4){
        convert(mat,produ,d,p,k);
        copy(mat,produ,n,n);

    }
    /*printf("INVERSE MATRIX:\n");
    print_mat(Ed,n);
    printf("\nINITIAL MATRIX:\n");
    print_mat(mat,n);
    */
    prod_ab(mat,n,n,Ed,n,n,produ);
    for(int i=0;i<n;i++){
        produ[i*n+i]-=1.0;
    }
    /*printf("\nRESULT:\n");
    print_mat(produ,n);
    */double tt=norm(produ,n);
    prod_ab(Ed,n,n,mat,n,n,produ);
    double tt1=norm(produ,n);
    if(tt>tt1)
        tt=tt1;
    printf("Residual: %e \n",tt);
    //sleep(2);
    delete[]produ;
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

double* read(FILE *f,int n){
    double *mat=new double [n*n];
    if(mat==NULL){
        printf("ERROR, no memory \n");
        return NULL;
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(fscanf(f,"%lf",&mat[i*n+j])!=1){
                printf("ERROR of fscanf\n");
                delete []mat;
                return NULL;
            }
        }
    }
    return mat;
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


double* read(FILE *f,int n,int d){

   // printf("%d azaza \n",*n);
    double temp;
    double *mat = new double[(n)*(n)];
    if(mat == NULL){                                        //very big matrix
        printf("NO MEMORY FOR MATRIX\n");

        return NULL;
    }
    int n_bl=(n) / d;
    int p=n%d;
   // printf("%d \n\n",n_bl);
    for(int i = 1; i <= (n); i++){
        for(int j = 1; j <= (n); j++){
            if(fscanf(f,"%lf",&temp)!=1){
                printf("ERROR, CAN'T READ MATRIX ELEMENTS\n");
                delete []mat;

                return NULL;
            }

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
    return mat;

}
void print_mat(double *a,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            printf("%10.10f ",a[i*n+j]);
        }
        printf("\n");

    }

}

