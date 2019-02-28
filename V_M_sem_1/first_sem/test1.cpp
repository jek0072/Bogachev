#include<stdio.h>
#include<math.h>

int main(){
    FILE *f;
    f=fopen("testik","w");
    if(f==NULL)
        return -1;
    int n=1000;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            fprintf(f,"%f ",fabs((double)(i-j)));
            printf("%f ",fabs((double)(i-j)));
        }
        fprintf(f,"\n");
        printf("\n");

    }
    return 0;
}
