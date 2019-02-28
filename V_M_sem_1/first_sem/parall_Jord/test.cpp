#include<stdio.h>
#include<time.h>
#include <stdlib.h>

int main(){
    FILE *f;
    f=fopen("in","w");
    
    srand(time(NULL));
    for(int i=0;i<3100;i++){
        for(int j=0;j<3100;j++){
            fprintf(f,"%d ",rand()%300);
        }
        fprintf(f,"\n");
    }
    return 0;
}
