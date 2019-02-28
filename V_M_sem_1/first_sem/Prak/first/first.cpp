#include<stdio.h>
#include<string.h>
#include<pthread.h>

#include<math.h>
#define EPS 1e-10


struct inform{
    char *file_name;
    pthread_mutex_t err_mut;

    pthread_mutex_t *mut;
    int* err;
    pthread_barrier_t *bar;
    int id;
    int max_id;
    double *max;
    int *max_ident;
    int* qq;
};











void* FUNC(void* T){


    inform *inf=(inform*)T;
    printf("%d\n",inf->id);
    FILE *f;
    f=fopen(inf->file_name,"r");
    if(f==NULL){
        pthread_mutex_lock(&(inf->err_mut));
        *inf->err=-1;
        pthread_mutex_unlock(&(inf->err_mut));
    }
    pthread_barrier_wait(inf->bar);

    pthread_mutex_lock(&(inf->err_mut));
        if(*inf->err<0){
            if(f!=NULL)
                fclose(f);

            pthread_mutex_unlock(&(inf->err_mut));
            return NULL;
        }
    pthread_mutex_unlock(&(inf->err_mut));
    int n=0;

    //rewind(f);
    double f1,f2;
    bool q=0;

    if(fscanf(f,"%lf",&f1)==1)
        n++;
    while(fscanf(f,"%lf",&f2)==1){
        n++;
        if(fabs(f1-f2)<EPS){
           if(q==0){
               inf->max[inf->id]=f1;

               inf->max_ident[inf->id]=1;
               q=1;
           }
           else{
               if(inf->max[inf->id]<f1)
                   inf->max[inf->id]=f1;
                   inf->max_ident[inf->id]=1;

           }
        }
        f1=f2;
    }





    pthread_barrier_wait(inf->bar);

    pthread_mutex_lock(&(inf->mut));
        if(inf->max_ident[inf->id]==1){
            if(inf->max_ident[0]==0){
                inf->max[0]=inf->max[inf->id];
                inf->max_ident[0]=1;
            }
            else{
                if(inf->max[0]<inf->max[inf->id])
                    inf->max[0]=inf->max[inf->id];
            }
        }
    pthread_mutex_unlock(&(inf->mut));

    pthread_barrier_wait(inf->bar);
    double max=inf->max[0];
   // printf("%f ",max);
    pthread_barrier_wait(inf->bar);
    rewind(f);

    int count=0;
    for(int i=0;i<n;i++){
        fscanf(f,"%lf",&f1);
        if(f1>max)
            count++;
    }



    pthread_mutex_lock(&(inf->mut));
   // if(inf->id==0)
    //    printf("%d HAHAHA",count);
    if(*inf->qq==0){
        *inf->qq=1;
        inf->max[0]=count;
    }
    else
        inf->max[0]+=count;
    pthread_mutex_unlock(&(inf->mut));

    fclose(f);
    //delete[]arr;
    return NULL;
}

















int main(int argc,char** argv){
    if(argc<=1){
        printf("USAGE: a.txt b.txt c.txt etc\n");
        return -1;
    }
    int n=argc-1;

    pthread_t* id=new pthread_t[n];
    if(id==NULL){
        printf("ERROR MEMORY\n");
        return -2;
    }

    inform *T=new inform[n];
    if(T==NULL){
        delete[] id;
        printf("ERROR memory\n");
        return -3;
    }

    int *err=new int[1];
    if(err==NULL){
        delete[]id;
        delete[]T;
        printf("ERROR memory\n");
        return -4;
    }
    double *max=new double[n];
    if(max==NULL){
        printf("ERROR memory\n");
        delete[]id;
        delete[]T;
        delete[]err;
        return -5;
    }
    int *max_ident=new int[n];
    if(max_ident==NULL){
        printf("ERROR memory\n");
        delete[]max;
        delete[]id;
        delete[]T;
        delete[]err;
        return -5;
    }
    int *qq=new int[1];
    for(int i=0;i<n;i++)
        max_ident[i]=0;

    pthread_mutex_t err_mut=PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t mut=PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier,NULL,n);
    *(err)=0;
    *(qq)=0;
    for(int i=0;i<n;i++){
        T[i].bar=&barrier;
        T[i].file_name=new char[100];
        strcpy(T[i].file_name,argv[i+1]);
        T[i].err_mut=err_mut;
        T[i].mut=&mut;
        T[i].err=err;
        T[i].id=i;
        T[i].qq=qq;
        T[i].max_id=n;
        T[i].max=max;
        T[i].max_ident=max_ident;
    }

    for(int i=0;i<n;i++){
        pthread_create(id+i,NULL,&FUNC,T+i);
    }

    for(int i=0;i<n;i++){
        pthread_join(id[i],NULL);
    }
    if(*err==-1){
      printf("CAN't open file\n");
      for(int i=0;i<n;i++)
          delete[]T[i].file_name;
      pthread_barrier_destroy(&barrier);
      delete[]err;
      delete[]max;
      delete[]max_ident;
      delete[]id;
      delete[]T;
      delete[]qq;
      return -11;
    }
    if(*err==-2){
        printf("NO constant sector\n");
        for(int i=0;i<n;i++)
            delete[]T[i].file_name;
        pthread_barrier_destroy(&barrier);
        delete[]err;
        delete[]max;
        delete[]max_ident;
        delete[]id;
        delete[]qq;
        delete[]T;
    }
    printf("HELLO\n");
   printf("%d\n",(int)max[0]);

    for(int i=0;i<n;i++)
        delete[]T[i].file_name;
    pthread_barrier_destroy(&barrier);
    delete[]err;
    delete[]max;
    delete[]max_ident;
    delete[]id;
    delete[]T;
    delete[]qq;
    return 0;
}
