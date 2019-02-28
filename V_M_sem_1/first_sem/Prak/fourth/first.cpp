#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<string.h>
#include<math.h>
# define EPS 1e-10
struct inform{
    double *a;
    int n;
    int t;
    int id;
    int beg;
    int end;
    int* kolvo;
    double *res;
    int *err;
    pthread_mutex_t *mutex;
    pthread_barrier_t *barrier;

};

void* func(void* T){
    inform *inf=(inform*)T;
    int n=inf->n;
    int t=inf->t;
    n=n;
    t=t;
    int beg=inf->beg;
    int end=inf->end;
    int id=inf->id;
    double *a=inf->a;
    int *kolvo=inf->kolvo;
    double *res=inf->res;

    n=n;
    t=t;
    kolvo=kolvo;
    res=res;
    beg=beg;
    id=id;
    end=end;
    a=a;
    for(int i=beg;i<end;i++)
    {
        if(i - 1 >= 0 && i + 1 < n)
        {
            if(a[i] <= a[i + 1] && a[i]  <= a[i - 1])
            {
                pthread_mutex_lock( inf->mutex);
                res[*kolvo]=a[i];
                (*kolvo)++;
                pthread_mutex_unlock( inf->mutex);
            }
        }
    }
    pthread_barrier_wait(inf->barrier);
    if(id==0){
        double sum=0.0;
        int kol=0;
        for(int i=0;i < *kolvo;i++){
            sum+=res[i];
            kol++;
        }
        if(kol == 0){
            *inf->err=-1;
            pthread_barrier_wait (inf->barrier);

            return NULL;
        }
        printf("ans = %f\n",sum/kol);
        res[0]=sum/kol;
        /*for(int i=0;i<*kolvo;i++){
            printf("%f ", a[i]);
        }*/
    }
    pthread_barrier_wait(inf->barrier);
    if(*inf->err == -1){
        return NULL;
    }
    int ri=0;

    int q1=0,q2=0;
    if(end-beg==1){
        ri=1;
        if(beg-1>=0 && beg+1 <n){
            if(a[beg]<a[beg-1] && a[beg]< a[beg+1]){
                q1=1;
            }
        }
    }

    if(end-beg>=2){
        ri=2;
        if(beg-1>=0 && beg+1 <n){
            if(a[beg]<a[beg-1] && a[beg]< a[beg+1]){
                q1=1;
            }
        }
        if(end-1 -1 >= 0 && end - 1 + 1 < n){
            if(a[end - 1]<a[end - 2] && a[end-1]< a[end]){
                q2=1;

            }
        }
    }
    pthread_barrier_wait(inf->barrier);

    bool fir=0,sec=0;
    for(int i=beg+1;i<end-1;i++){
        if(a[i]<=a[i+1] && a[i] <= a[i-1]){
            sec=1;
        }
        else
            sec=0;
        if(fir == 1){
            a[i-1]=res[0];
        }
        fir=sec;
    }
    if(fir == 1)
        a[end-2]=res[0];

    pthread_barrier_wait(inf->barrier);
    if(ri==1){
        if(q1==1)
            a[beg]=res[0];
    }
    if(ri==2){
        if(q1==1)
            a[beg]=res[0];
        if(q2==1)
            a[end-1]=res[0];
    }
    pthread_barrier_wait(inf->barrier);

    if(id==0){
        printf("res: ");
        for(int i=0;i<n;i++){
                    printf("%f ", a[i]);
        }
        printf("\n");
    }

    return NULL;

}
int main(int argc, char **argv){
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
    if(t<0){
        return -3123;
    }
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
            if(f!=NULL)
                fclose(f);
            return -2;

        }
    }
    pthread_t *id=new pthread_t[t];
    if(id==NULL){
        printf("ERROR id\n");
        delete[] ar;
        if(f!=NULL)
            fclose(f);
        return -4;

    }

    inform *inf=new inform[t];
    if(inf==NULL){
        printf("ERROR inf\n");
        delete[] ar;
        delete[] id;
        if(f!=NULL)
            fclose(f);
        return -5;

    }

    pthread_mutex_t mut=PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier,NULL,t);
    int kol=0;
    int err=0;
    double *res=new double [n];
    for(int i=0;i<t-1;i++){
           inf[i].a=ar;
           inf[i].n=n;
           inf[i].id=i;
           inf[i].t=t;
           inf[i].err=&err;
           inf[i].beg=(n/t)*i;
           inf[i].end=(n/t)*i+n/t;
           inf[i].barrier=&barrier;
           inf[i].mutex=&mut;
           inf[i].res=res;
           inf[i].kolvo=&kol;

    }
    inf[t-1].beg=(n/t)*(t-1);
    inf[t-1].end=n;
    inf[t-1].err=&err;

    inf[t-1].a=ar;
    inf[t-1].n=n;
    inf[t-1].id=t-1;
    inf[t-1].t=t;
    inf[t-1].barrier=&barrier;
    inf[t-1].mutex=&mut;
    inf[t-1].res=res;
    inf[t-1].kolvo=&kol;

    int tttt;
    for(int i=0;i<t;i++){
        tttt=pthread_create(id+i,NULL,&func,inf+i);
        if(tttt<0)
        {
            printf("ERROR of create\n");
            if(f!=NULL)
                fclose(f);
            pthread_barrier_destroy(&barrier);
            delete[]res;
            delete[]inf;
            delete[]ar;
            delete[]id;
            return -6;
        }
    }
    for(int i=0;i<t;i++)
        pthread_join(id[i],NULL);
    if(err==-1){
        printf("no local min\n");
    }
    if(f!=NULL)
        fclose(f);
    pthread_barrier_destroy(&barrier);
    delete[]res;
    delete[]inf;
    delete[]ar;
    delete[]id;
    return 0;

}


