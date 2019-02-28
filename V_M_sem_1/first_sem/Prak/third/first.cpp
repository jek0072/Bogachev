#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<pthread.h>

#include <sys/time.h>
#include <sys/resource.h>

struct inform{
    double *ar;
    int n;
    int p;
    int id;
    int beg;
    int end;
    double *kray;
    double *kray_ident;
    pthread_barrier_t *bar;
    pthread_mutex_t *mut;
};

double get_time(){
    struct rusage buf;
    getrusage(RUSAGE_THREAD,&buf);
    return (double)buf.ru_utime.tv_sec+(double)buf.ru_utime.tv_usec/1000000.;
}

double get_full_time(){
    struct timeval buf;
    gettimeofday(&buf,0);
    return (double)buf.tv_sec+(double)buf.tv_usec/1000000.;
}

void* func(void* T){
    inform* inf=(inform*)T;
    int n=inf->n;
    int p=inf->p;
    n=n;
    p=p;
    int beg=inf->beg;
    int end=inf->end;
    int id=inf->id;
    double *ar=inf->ar;
    beg=beg;
    id=id;
    end=end;
    ar=ar;

    /*pthread_mutex_lock(inf->mut);
        printf("id: %d\n",id);
        for(int i=beg;i<end;i++){
            printf("%f ",ar[i]);
        }
        printf("\n");
    pthread_mutex_unlock(inf->mut);
*/


    int count=end-beg;
    double temp[4];
    int temp_1[4];
    for(int i=0;i<4;i++){
        temp[i]=-20;
        temp_1[i]=-20;
    }

    double t1=get_time();
    int ident=-1;
    if(count==1){

        if(beg-2>=0 && beg+2<n){
            temp[0]=(ar[beg-2]+ar[beg+2])/2.0;
            temp_1[0]=1;
            ident=1;

        }
    }

    if(count==2){

        if(beg-2>=0 && beg+2<n){
            temp[0]=(ar[beg-2]+ar[beg+2])/2.0;
            temp_1[0]=1;
            ident=2;

        }

        if(beg+1-2>=0 && beg+1+2<n){
            temp[1]=(ar[beg-2+1]+ar[beg+2+1])/2.0;
            temp_1[1]=1;
            ident=2;

        }

    }
    if(count==3){
        if(beg-2>=0 && beg+2<n){
            temp[0]=(ar[beg-2]+ar[beg+2])/2.0;
            temp_1[0]=1;
            ident=3;

        }

        if(beg+1-2>=0 && beg+1+2<n){
            temp[1]=(ar[beg-2+1]+ar[beg+2+1])/2.0;
            temp_1[1]=1;
            ident=3;

        }
        if(beg+2-2>=0 && beg+2+2<n){
            temp[2]=(ar[beg-2+2]+ar[beg+2+2])/2.0;
            temp_1[2]=1;
            ident=3;

        }

    }

    if(count>=4){
        if(beg-2>=0 && beg+2<n){
            temp[0]=(ar[beg-2]+ar[beg+2])/2.0;
            temp_1[0]=1;
            ident=4;

        }

        if(beg+1-2>=0 && beg+1+2<n){
            temp[1]=(ar[beg-2+1]+ar[beg+2+1])/2.0;
            temp_1[1]=1;
            ident=4;

        }
        if(end-2-2>=0 && end-2+2<n){
            temp[2]=(ar[end-2-2]+ar[end-2+2])/2.0;
            temp_1[2]=1;
            ident=4;

        }

        if(end-2-1>=0 && end+2-1<n){
            temp[3]=(ar[end-1-2]+ar[end-1+2])/2.0;
            temp_1[3]=1;
            ident=4;

        }
    }
    /*pthread_mutex_lock(inf->mut);
        printf("id=%d\n %f %f %f %f\n",id,temp[0],temp[1],temp[2],temp[3]);
    pthread_mutex_unlock(inf->mut);
*/
    pthread_barrier_wait(inf->bar);
    int help[3];
    help[0]=0;
    help[1]=0;
    help[2]=0;

    double help_1[3];
    help_1[0]=0;
    help_1[1]=0;
    help_1[2]=0;

    for(int i=beg+2;i<end-2;i++){


        if(i-2 >= 0 && i+2 < n){
            help[2]=1;
            help_1[2]=(ar[i-2]+ar[i+2])/2.0;

        }
        if(help[0]==1){
            ar[i-2]=help_1[0];
        }

        help[0]=help[1];
        help[1]=help[2];
        help_1[0]=help_1[1];
        help_1[1]=help_1[2];

    }
    if(help[0]==1){
        ar[end-4]=help_1[0];
    }
    if(help[1]==1){
        ar[end-3]=help_1[1];
    }

    pthread_barrier_wait(inf->bar);
    switch(ident){
        case 1:{
            ar[beg]=temp[0];
            break;
        }
        case 2:{
            if(temp_1[0]==1){
                ar[beg]=temp[0];
            }
            if(temp_1[1]==1){
                ar[beg+1]=temp[1];
            }
            break;
        }
        case 3:{
            if(temp_1[0]==1){
                ar[beg]=temp[0];
            }
            if(temp_1[1]==1){
                ar[beg+1]=temp[1];
            }
            if(temp_1[2]==1){
                ar[beg+2]=temp[2];
            }
            break;
        }

        case 4:{
            if(temp_1[0]==1){
                ar[beg]=temp[0];
            }
            if(temp_1[1]==1){
                ar[beg+1]=temp[1];
            }
            if(temp_1[2]==1){
                ar[end-2]=temp[2];
            }
            if(temp_1[3]==1){
                ar[end-1]=temp[3];
            }
        }


    }
    //pthread_barrier_wait(inf->bar);

    t1=get_time()-t1;

    pthread_mutex_lock(inf->mut);
        printf("id=%d: %f\n",id,t1);
    pthread_mutex_unlock(inf->mut);
    return NULL;
}


int main(int argc, char ** argv){
    int n;
    int t=0;
    char file_name[100];
    if(argc!=4 ){
        printf("USAGE: ./a.out n t input.txt\n");
        return -1;

    }
    strcpy(file_name,argv[3]);
    n=atoi(argv[1]);
    t=atoi(argv[2]);
    if(t>n)
        t=n-1;

    //printf("%d %d %s\n",n,t,file_name);
    FILE *f;
    f=fopen(file_name,"r");
    if(f==NULL){
        printf("FILE not exist\n");
        return -2;
    }
    double *ar=new double[n];
    if(ar==NULL){
        printf("ERROR of memory\n");
        return -3;

    }

    for(int i=0;i<n;i++){
        if(fscanf(f,"%lf",&ar[i])!=1){
            printf("ERROR read\n");
            delete[] ar;
            return -2;

        }
    }
    pthread_t *id=new pthread_t[t];
    if(id==NULL){
        printf("ERROR id\n");
        delete[] ar;
        return -4;

    }

    inform *inf=new inform[t];
    if(inf==NULL){
        printf("ERROR inf\n");
        delete[] ar;
        delete[] id;
        return -5;

    }

    pthread_mutex_t mut=PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier,NULL,t);

    for(int i=0;i<t-1;i++){
        inf[i].id=i;
        inf[i].n=n;
        inf[i].p=t;
        inf[i].ar=ar;
        inf[i].beg=(n/t)*i;
        inf[i].end=(n/t)*i+n/t;
        inf[i].mut=&mut;
        inf[i].bar=&barrier;
    }
    inf[t-1].id=t-1;
    inf[t-1].n=n;
    inf[t-1].p=t;
    inf[t-1].ar=ar;
    inf[t-1].beg=(n/t)*(t-1);
    inf[t-1].end=n;
    inf[t-1].mut=&mut;
    inf[t-1].bar=&barrier;
    int temp;
    /*for(int i=0;i<n;i++){
        printf("%f  ",ar[i]);

    }*/
    printf("\n");
    double t0=get_full_time();
    for(int i=0;i<t;i++){
        temp=pthread_create(id+i,NULL,&func,inf+i);
        if(temp<0){

            pthread_barrier_destroy(&barrier);
            delete[] ar;
            delete[]id;
            delete[] inf;
            return -6;
        }
    }
    for(int i=0;i<t;i++){
        pthread_join(id[i],NULL);

    }

    t0=get_full_time()-t0;
    printf("ELAPSED: %f\n",t0);

    int n1;
    if(n>30)
        n1=30;
    else
        n1=n;
    printf ("Result: ");
    for(int i=0;i<n1;i++){
        printf(" %.2f",ar[i]);

    }
    printf("\n");

    pthread_barrier_destroy(&barrier);
    delete[] ar;
    delete[]id;
    delete[] inf;
    return 0;
}
