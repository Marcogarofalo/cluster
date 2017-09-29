#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


double main(int argc, char* argv[]){
 double ave=0,sigma=0;
 int i;
 
 for(i=1;i<argc;i++)
     ave+=atof(argv[i]);
 ave/=(double)(argc-1);
 for(i=1;i<argc;i++)
     sigma+=(ave-atof(argv[i]))*(ave-atof(argv[i]));
/* printf("%f\n",atof(argv[1]));*/
 sigma/=(double)(argc-2);
 sigma=sqrt(sigma/((double)(argc-1)));
 printf("%f   %f",ave,sigma);

 return ave;
}
