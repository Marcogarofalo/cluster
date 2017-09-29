#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


double main(int argc, char* argv[]){
 double lambda,kappa;
 int i;

if(argc!=3) {printf("usage: ./lambda $kappa $lambda \n"), exit(0);}
 
kappa=atof(argv[1]);
lambda=atof(argv[2]);
lambda=4*kappa*kappa*lambda;
printf("%.6f",lambda);

 return lambda;
}
