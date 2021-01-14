#include <array>
#include <cmath>
#include <ctime>
#include <vector>

#include "mdp.h"


int main(int argc, char** argv){

    int T, L;
    if (argc!=5){
       printf("usage:  ./generate_trivial_scalar  -T T  -L  L \n");
    }
    for(int i=0; i<argc; i++  ){
       if (strcmp(argv[i],"-T")==0 )
          T=atoi(argv[i+1]);
       if (strcmp(argv[i],"-L")==0 )
          L=atoi(argv[i+1]);
    }
    printf("T = %d   L=%d \n",T,L);

    double *phi=(double*) calloc(T*L*L*L*4, sizeof(double));
   
    FILE *f0=NULL;
    f0=fopen("zero_scalar","w+");
    if (f0==NULL) { printf("error opening zero_scalar\n"); exit(1);}
    fwrite(phi,sizeof(double), T*L*L*L*4, f0 );
    fclose(f0);
 
    FILE *f1=NULL;
    f1=fopen("identity_scalar","w+");
    for (size_t x=0; x<T*L*L*L ; x++)
        phi[x*4]=1.;

    if (f1==NULL) { printf("error opening identity_scalar\n"); exit(1);}
    fwrite(phi,sizeof(double), T*L*L*L*4, f1 );
    fclose(f1);
     
  
}
