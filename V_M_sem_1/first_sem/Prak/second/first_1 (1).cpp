#include<stdio.h>
#include<string.h>
#include<pthread.h>

#include<math.h>
#define EPS 1e-10


struct inform{
    char *file_name;
    pthread_mutex_t* err_mut;
	int *res;
    pthread_mutex_t* mut;
    int* err;
    pthread_barrier_t *bar;
    int id;
    int max_id;
    double *share;
    double *sum_el;
    int *num_el;
    int *share_help;
};


int determ_type(int *a){
    if(a[0]==0 && a[1]==0 && a[2]==0 && a[3]==0)
        return 0;
    if(a[0]==1 && a[1]==1 && a[2]==0 && a[3]==0)
        return 1;
    if(a[0]==1 && a[1]==1 && a[2]==1 && a[3]==1)
        return 2;
    return -1;
}








void* FUNC(void* T){


    inform *inf=(inform*)T;
   // printf("%d\n",inf->id);
    FILE *f;
    f=fopen(inf->file_name,"r");
    if(f==NULL){
        pthread_mutex_lock((inf->err_mut));
        *inf->err=-1;
        pthread_mutex_unlock((inf->err_mut));
    }
    pthread_barrier_wait(inf->bar);

    pthread_mutex_lock((inf->err_mut));
        if(*inf->err<0){
            if(f!=NULL)
                fclose(f);

            pthread_mutex_unlock((inf->err_mut));
            return NULL;
        }
    pthread_mutex_unlock((inf->err_mut));


    double f1=0,f2=0,f3=0;//                                              scanf three numbers and check them
    int t_tt=fscanf(f,"%lf %lf %lf",&f1,&f2,&f3);

    if(t_tt<=0){
        inf->share_help[4*inf->id]=0;
        inf->share_help[4*inf->id+1]=0;

        inf->share_help[4*inf->id+2]=0;
        inf->share_help[4*inf->id+3]=0;
    }

    if(t_tt==1){
        inf->share_help[4*inf->id]=1;
        inf->share_help[4*inf->id+1]=1;
        inf->share_help[4*inf->id+2]=0;
        inf->share_help[4*inf->id+3]=0;
        inf->share[4*inf->id]=f1;
        inf->share[4*inf->id+1]=f1;
		inf->share[4*inf->id+2]=f1;
		inf->share[4*inf->id+3]=f1;

    }

    if(t_tt==2){
        inf->share_help[4*inf->id]=1;
        inf->share_help[4*inf->id+1]=1;
        inf->share_help[4*inf->id+2]=1;
        inf->share_help[4*inf->id+3]=1;
        inf->share[4*inf->id]=f1;
        inf->share[4*inf->id+1]=f2;
        inf->share[4*inf->id+2]=f1;
        inf->share[4*inf->id+3]=f2;
    }

    if(t_tt==3){
        inf->share_help[4*inf->id]=1;
        inf->share_help[4*inf->id+1]=1;
        inf->share_help[4*inf->id+2]=1;
        inf->share_help[4*inf->id+3]=1;
        inf->share[4*inf->id]=f1;
        inf->share[4*inf->id+1]=f2;
    }

    double sum=0;
    int amount=0;
    if(t_tt==3){

       /////////////////////
        do{
            if(f2<=f1 && f2<=f3){
                sum+=f2;
                amount++;
            }
            f1=f2;
            f2=f3;
        }while(fscanf(f,"%lf",&f3)==1);
        ///////////////////////////


            //printf("HHHHHHH %f %f\n",f1,f2);
            f3=f2;
            f2=f1;
            //printf("HHHHHHH %f %f\n",f2,f3);
            inf->share[4*inf->id+2]=f2;
            inf->share[4*inf->id+3]=f3;

    }
    pthread_barrier_wait(inf->bar);//barier for full search
	
	
    pthread_mutex_lock(inf->mut);//////////////////////////////
    *inf->num_el+=amount;
    *inf->sum_el+=sum;
    //printf("AAAAAAAAZAAAAAAAA %f %d\n",sum,amount);
    //printf("Id=%d     %d %f\n",inf->id, amount, sum);
    pthread_mutex_unlock(inf->mut);////////////////////////
	
	if(inf->id==0){
		
		int n_new=0;
        amount=0;
        sum=0;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int j=0;
		for(int i=0;i<inf->max_id;i++){
			if(determ_type(inf->share_help+4*i)==0){
				
				if(i>=j)
					j=i+1;
				while(j<inf->max_id && determ_type(inf->share_help+4*j)==0)
				{	j++;
					
				}
				if(j < inf->max_id){
						
					inf->share[4*i]=inf->share[4*j];
					inf->share[4*i+1]=inf->share[4*j+1];
					inf->share[4*i+2]=inf->share[4*j+2];
					inf->share[4*i+3]=inf->share[4*j+3];
					inf->share_help[4*i]=inf->share_help[4*j];
					inf->share_help[4*i+1]=inf->share_help[4*j+1];
					inf->share_help[4*i+2]=inf->share_help[4*j+2];
					inf->share_help[4*i+3]=inf->share_help[4*j+3];
					inf->share_help[4*j]=0;
					inf->share_help[4*j+1]=0;
					inf->share_help[4*j+2]=0;
					inf->share_help[4*j+3]=0;
				}
				else
					break;
			}
				n_new++;
		}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for(int i=1;i<n_new-1;i++){
			if(determ_type(inf->share_help+i*4)==2){
				
					if(inf->share[i*4]<=inf->share[i*4+1] && inf->share[i*4]<=inf->share[i*4-1]){//left
						amount++;
						sum+=inf->share[4*i];
						
					}
					
					
					
					
					///////////////////////////////////////////////////////////////////////////////////////////////////
					if(inf->share[i*4+3]<=inf->share[i*4+2] && inf->share[i*4+3]<=inf->share[i*4+4]){//right
						amount++;
						sum+=inf->share[4*i+3];
						
					}
					
			}
			else{//type 1
				if(inf->share[i*4]<=inf->share[i*4+4] && inf->share[i*4]<=inf->share[i*4-1]){
					amount++;
					sum+=inf->share[4*i];
					
				}
			}
		}
		
		if(n_new>=2){
			if(determ_type(inf->share_help)==2){
				if(inf->share[3]<=inf->share[2] && inf->share[3]<=inf->share[4]){//right
						amount++;
						sum+=inf->share[3];
					
					}
			}
			if(determ_type(inf->share_help+(n_new-1)*4)==2){
				if(inf->share[(n_new-1)*4]<=inf->share[(n_new-1)*4+1] && inf->share[(n_new-1)*4]<=inf->share[(n_new-1)*4-1]){//left
					amount++;
					sum+=inf->share[(n_new-1)*4];
				
				}
			}
		}
		pthread_mutex_lock(inf->mut);//////////////////////////////
		*inf->num_el+=amount;
		*inf->sum_el+=sum;
		//printf("Id=%d     %d %f\n",inf->id, amount, sum);
		pthread_mutex_unlock(inf->mut);////////////////////////
		
    }
    
	pthread_barrier_wait(inf->bar);
	rewind(f);
	int result=0;
	while(fscanf(f,"%lf",&f1)==1){
		if(f1<(*inf->sum_el)/(*inf->num_el))
			result++;
	//	printf("%d %d\n",inf->id,result);
	}
	
	pthread_mutex_lock(inf->mut);//////////////////////////////
	*inf->res+=result;
	pthread_mutex_unlock(inf->mut);////////////////////////
	
	fclose(f);
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


    double *share=new double[4*n];
    if(share==NULL){
        printf("ERROR memory\n");
        delete[]id;
        delete[]T;
        delete[]err;
        return -5;
    }


    int *share_help=new int[4*n];
    if(share_help==NULL){
        printf("ERROR memory\n");
        delete[]id;
        delete[]T;
        delete[]err;
        delete[] share;
        return -5;
    }
	int *res=new int[1];
	if(res==NULL){
		
        printf("ERROR memory\n");
        delete[]id;
        delete[]T;
        delete[]err;
        delete[] share;
        
        delete[] share_help;
        return -6;
	}
	
    for(int i=0;i<4*n;i++){
        share_help[i]=0;
    }
    double *sum_el=new double[1];
    int* num_el=new int[1];
    *sum_el=0;
    *num_el=0;
    *res=0;
    pthread_mutex_t err_mut=PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t mut=PTHREAD_MUTEX_INITIALIZER;
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier,NULL,n);
    *(err)=0;
    for(int i=0;i<n;i++){
        T[i].bar=&barrier;
        T[i].file_name=new char[100];
        strcpy(T[i].file_name,argv[i+1]);
        T[i].err_mut=&err_mut;
        T[i].mut=&mut;
        T[i].sum_el=sum_el;
        T[i].num_el=num_el;
        T[i].err=err;
        T[i].id=i;
        T[i].max_id=n;
        T[i].share=share;
        T[i].share_help=share_help;
		T[i].res=res;
    }
    //printf("  %d  ",n);
    for(int i=0;i<n;i++){
        pthread_create(id+i,NULL,&FUNC,T+i);
    }

    for(int i=0;i<n;i++){
        pthread_join(id[i],NULL);
    }

    if(*err==-1){
      printf("HERE CAN't open file\n");

      for(int i=0;i<n;i++)
            delete[]T[i].file_name;
      pthread_barrier_destroy(&barrier);
      delete[]err;
      delete[]share;

      delete[] sum_el;
	  delete[] res;
      delete[] num_el;
      delete[] share_help;
      delete[]id;
      delete[]T;
      return -11;
    }
    printf("RES: %f %d\n",*sum_el,*num_el);
	printf("Ans: %d\n",*res);

     delete[] share_help;
	
    for(int i=0;i<n;i++)
        if(T[i].file_name!=NULL)
            delete[]T[i].file_name;
    pthread_barrier_destroy(&barrier);
    delete[]err;
    delete[]share;
    delete[] sum_el;
	
	delete[] res;
    delete[] num_el;
    delete[]id;
    delete[]T;
    return 0;


}

