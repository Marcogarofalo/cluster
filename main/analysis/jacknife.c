#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

//#include <cmath>
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <iostream>
//#include <getopt.h>

//using namespace std;

// global lattice data
int L0, L1, L2, L3, V,global_c=1;
int verbose;
char* input_filename, *input_filename1;
char* output_filename;
int fit_min=2,fit_max=10;

double* phi_field;

double **gammaFbb, *CbbF,*taubb_intF,*dtau,*dobs,*ddobs,*obs,*abb;
double ***corr_ave, ***effective_mass, ***matrix_element;
int *w,N_tot,every,call_f,global_k=0;
//static const std::complex<double> I (0.0, 1.0);
double *constant_fit(int Lmin, int Lmax ,double **s)
{   
    int i;
    double *r;
    
    r=(double*) calloc(2,sizeof(double));
    
    for(i=Lmin;i<Lmax;i++){
        r[0]+=s[0][i]/(s[1][i]*s[1][i]);
        r[1]+=1./(s[1][i]*s[1][i]);
    }
    r[0]/=r[1];
    r[1]=sqrt(1./r[1]);
    
    return r;
    
}


double **Gamma(int t, int var, int order, int rep,int nconf, double *a, double *bba)
{
  double **r;
  int i0,i1,i2,i3,alpha,N;
  
  alpha=(order+1)*var;
  N=alpha*rep;
  r=(double**) malloc(sizeof(double*)*(alpha));
  for(i0=0;i0<alpha;i0++)
    r[i0]=(double*) calloc(alpha,sizeof(double));
  
  for(i0=0;i0<alpha;i0++)
    for(i1=0;i1<alpha;i1++)
        for(i2=0;i2<rep;i2++)
            for(i3=0;i3<(nconf-t);i3++)
                r[i0][i1]+= (a[i0+i2*alpha+i3*N]-bba[i0])*(a[i1+i2*alpha+(i3+t)*N]-bba[i1]);
    // if(t==0) printf("r[0][0]=%f\n",r[0][0]);  
          
    
          /*printf("a-bba=%f\n",a[0]-bba[0]);*/
for(i0=0;i0<alpha;i0++)
  for(i1=0;i1<alpha;i1++)
	    r[i0][i1]/=((double)(rep*nconf-rep*t));
     
return r;
}


double *function(int var, int order,int flow ,double *ah)
{
double *r,u,d,ct[2],ctp[2],ctm[2],mass,tmp_mass,res;
int i,k;
r=(double*) calloc((order+1),sizeof(double)); 
if(call_f==-1){
    for(i=10;i<19;i++){
       ct[0]=(ah[0+i*2]+ah[0+(L0-i)*2])/2.;
       ct[1]=ah[1];
       ct[0]=ct[0]-ct[1]*ct[1];
	
       ctp[0]=(ah[0+(i+1)*2]+ah[0+(L0-(i+1))*2])/2.;
       ctp[1]=(ah[1+(i+1)*2]+ah[1+(L0-(i+1))*2])/2.;
       ctp[0]=ctp[0]-ct[1]*ct[1];


        mass=log(ct[0]/ctp[0]);
        res=1;
        while(res>1e-15){
            u=1.+exp(-mass*(L0-2*i-2));
            d=1.+exp(-mass*(L0-2*i));
            tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
            res=tmp_mass-mass;
            mass=tmp_mass;
        }
        r[0]+=mass;
        
    }
    r[0]/=9.;
}
else if(call_f==1) {
     i=flow;
     if(i !=0){
        ct[0]=(ah[0+i*2]+ah[0+(L0-i)*2])/2.;
        ct[1]=(ah[1+i*2]+ah[1+(L0-i)*2])/2.;
        ct[0]=ct[0]-ah[1]*ah[1];
        r[0]=ct[0];
       
     } else{
        ct[0]=ah[0+i*2];
        ct[1]=ah[1+i*2];
        ct[0]=ct[0]-ah[1]*ah[1];
        r[0]=ct[0];  
    }
    
}
else if(call_f==3) {
    double **masses,*fit;
    
    masses=(double**) malloc(sizeof(double*)*2);
    masses[0]=(double*) malloc(sizeof(double)*L0);
    masses[1]=(double*) malloc(sizeof(double)*L0);
    
    
    for(i=0;i<L0/2;i++){
        if(i!=0){
            ct[0]=(ah[0+i*2]+ah[0+(L0-i)*2])/2.;
            ct[1]=(ah[1+i*2]+ah[1+(L0-i)*2])/2.;
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=(ah[0+(i+1)*2]+ah[0+(L0-(i+1))*2])/2.;
            ctp[1]=(ah[1+(i+1)*2]+ah[1+(L0-(i+1))*2])/2.;
            ctm[0]=(ah[0+(i-1)*2]+ah[0+(L0-(i-1))*2])/2.;
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];

        }
        else {
            ct[0]=ah[0+i*2];
            ct[1]=ah[1+i*2];
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=ah[0+(i+1)*2];
            ctp[1]=ah[1+(i+1)*2];
            ctm[0]=ah[0+(i-1)*2];
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];
        }
        
        mass=log(ct[0]/ctp[0]);

        res=1;
        while(res>1e-19){
            u=1.+exp(-mass*(L0-2*i-2));
            d=1.+exp(-mass*(L0-2*i));
            tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
            res=tmp_mass - mass;
            mass=tmp_mass;
        }
        
        
        masses[0][i]=mass;
        masses[1][i]=effective_mass[global_k][1][i];
    }
    
     fit=constant_fit(fit_min,fit_max,masses);
  
     
     ct[0]=(ah[0+flow*2]+ah[0+(L0-flow)*2])/2.;
     if(flow==0)ct[0]=ah[0+flow*2];
     
     r[0]=ct[0]/(exp(-fit[0]*flow)+exp(-fit[0]*(L0-flow)));
 //  if (flow==0 && global_c==1){ printf("%f   %f  %.15g  \n",fit[0],ct[0],(exp(-fit[0]*flow)+exp(-fit[0]*(L0-flow)))); global_c=2;}
     free(masses[0]);free(masses[1]);free(fit);free(masses);
     
}

else if(call_f==2) {
     if(flow!=0){
        ct[0]=(ah[0+flow*2]+ah[0+(L0-flow)*2])/2.;
        ct[1]=(ah[1+flow*2]+ah[1+(L0-flow)*2])/2.;
            ct[0]=ct[0]-ah[1]*ah[1];

        ctp[0]=(ah[0+(flow+1)*2]+ah[0+(L0-(flow+1))*2])/2.;
        ctp[1]=(ah[1+(flow+1)*2]+ah[1+(L0-(flow+1))*2])/2.;
        ctm[0]=(ah[0+(flow-1)*2]+ah[0+(L0-(flow-1))*2])/2.;
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];

     }
     else {
        ct[0]=ah[0+flow*2];
        ct[1]=ah[1+flow*2];
            ct[0]=ct[0]-ah[1]*ah[1];

        ctp[0]=ah[0+(flow+1)*2];
        ctp[1]=ah[1+(flow+1)*2];
        ctm[0]=ah[0+(flow-1)*2];
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];
     }
    
     mass=log(ct[0]/ctp[0]);

     res=1;
     while(res>1e-16){
        u=1.+exp(-mass*(L0-2*flow-2));
        d=1.+exp(-mass*(L0-2*flow));
        tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
        res=tmp_mass - mass;
        mass=tmp_mass;
     }

     r[0]=mass;//acosh( (ctp[0]+ctm[0])/(2*ct[0])   );
}
 

else if(call_f==4) {
    double **masses,*fit;
    
    masses=(double**) malloc(sizeof(double*)*2);
    masses[0]=(double*) malloc(sizeof(double)*L0);
    masses[1]=(double*) malloc(sizeof(double)*L0);
    
    
    for(i=0;i<L0/2;i++){
        if(i!=0){
            ct[0]=(ah[0+i*2]+ah[0+(L0-i)*2])/2.;
            ct[1]=(ah[1+i*2]+ah[1+(L0-i)*2])/2.;
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=(ah[0+(i+1)*2]+ah[0+(L0-(i+1))*2])/2.;
            ctp[1]=(ah[1+(i+1)*2]+ah[1+(L0-(i+1))*2])/2.;
            ctm[0]=(ah[0+(i-1)*2]+ah[0+(L0-(i-1))*2])/2.;
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];

        }
        else {
            ct[0]=ah[0+i*2];
            ct[1]=ah[1+i*2];
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=ah[0+(i+1)*2];
            ctp[1]=ah[1+(i+1)*2];
            ctm[0]=ah[0+(i-1)*2];
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];
        }
        
        mass=log(ct[0]/ctp[0]);

        res=1;
        while(res>1e-19){
            u=1.+exp(-mass*(L0-2*i-2));
            d=1.+exp(-mass*(L0-2*i));
            tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
            res=tmp_mass - mass;
            mass=tmp_mass;
        }
        
        
        masses[0][i]=mass;
        masses[1][i]=effective_mass[global_k][1][i];
    }
    
     fit=constant_fit(fit_min,fit_max,masses);
  
     
//double p=2*sin(3.1415926535/L0 );
//p*=p;
     
    r[0]=fit[0]*fit[0]*( ah[1+1*2]-((double)V)* ah[1]*ah[1] ) ;
    //r[0]=ah[1+2*2]*(p+fit[0]*fit[0]);
 //  if (flow==0 && global_c==1){ printf("%f   %f  %.15g  \n",fit[0],ct[0],(exp(-fit[0]*flow)+exp(-fit[0]*(L0-flow)))); global_c=2;}
     free(masses[0]);free(masses[1]);free(fit);free(masses);
     
}

else if(call_f==7) {
    double **masses,*fit;
    
    masses=(double**) malloc(sizeof(double*)*2);
    masses[0]=(double*) malloc(sizeof(double)*L0);
    masses[1]=(double*) malloc(sizeof(double)*L0);
    
    
    for(i=0;i<L0/2;i++){
        if(i!=0){
            ct[0]=(ah[0+i*2]+ah[0+(L0-i)*2])/2.;
            ct[1]=(ah[1+i*2]+ah[1+(L0-i)*2])/2.;
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=(ah[0+(i+1)*2]+ah[0+(L0-(i+1))*2])/2.;
            ctp[1]=(ah[1+(i+1)*2]+ah[1+(L0-(i+1))*2])/2.;
            ctm[0]=(ah[0+(i-1)*2]+ah[0+(L0-(i-1))*2])/2.;
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];

        }
        else {
            ct[0]=ah[0+i*2];
            ct[1]=ah[1+i*2];
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=ah[0+(i+1)*2];
            ctp[1]=ah[1+(i+1)*2];
            ctm[0]=ah[0+(i-1)*2];
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];
        }
        
        mass=log(ct[0]/ctp[0]);

        res=1;
        while(res>1e-19){
            u=1.+exp(-mass*(L0-2*i-2));
            d=1.+exp(-mass*(L0-2*i));
            tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
            res=tmp_mass - mass;
            mass=tmp_mass;
        }
        
        
        masses[0][i]=mass;
        masses[1][i]=effective_mass[global_k][1][i];
    }
    
     fit=constant_fit(fit_min,fit_max,masses);
  
     
     
     r[0]=fit[0]*fit[0] ;
 //  if (flow==0 && global_c==1){ printf("%f   %f  %.15g  \n",fit[0],ct[0],(exp(-fit[0]*flow)+exp(-fit[0]*(L0-flow)))); global_c=2;}
     free(masses[0]);free(masses[1]);free(fit);free(masses);
     
}
else if(call_f==5) {
    double **masses,*fit;
    
    masses=(double**) malloc(sizeof(double*)*2);
    masses[0]=(double*) malloc(sizeof(double)*L0);
    masses[1]=(double*) malloc(sizeof(double)*L0);
    
    
    for(i=0;i<L0/2;i++){
        if(i!=0){
            ct[0]=(ah[0+i*2]+ah[0+(L0-i)*2])/2.;
            ct[1]=(ah[1+i*2]+ah[1+(L0-i)*2])/2.;
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=(ah[0+(i+1)*2]+ah[0+(L0-(i+1))*2])/2.;
            ctp[1]=(ah[1+(i+1)*2]+ah[1+(L0-(i+1))*2])/2.;
            ctm[0]=(ah[0+(i-1)*2]+ah[0+(L0-(i-1))*2])/2.;
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];

        }
        else {
            ct[0]=ah[0+i*2];
            ct[1]=ah[1+i*2];
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=ah[0+(i+1)*2];
            ctp[1]=ah[1+(i+1)*2];
            ctm[0]=ah[0+(i-1)*2];
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];
        }
        
        mass=log(ct[0]/ctp[0]);

        res=1;
        while(res>1e-19){
            u=1.+exp(-mass*(L0-2*i-2));
            d=1.+exp(-mass*(L0-2*i));
            tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
            res=tmp_mass - mass;
            mass=tmp_mass;
        }
        
        
        masses[0][i]=mass;
        masses[1][i]=effective_mass[global_k][1][i];
    }
    
     fit=constant_fit(fit_min,fit_max,masses);
  
     
//double p=2*sin(3.1415926535/L0 );
//p*=p;
     
    //r[0]=ah[1]*ah[1]/(ah[1+2*2]*(p+fit[0]*fit[0]));
     
     r[0]=ah[1]*ah[1]/(fit[0]*fit[0]*( ah[1+1*2]-((double)V)* ah[1]*ah[1] )) ;
//     r[0]=(ah[1+1*2]/(fit[0]*fit[0]*( ah[1+1*2]-((double)V)* ah[1]*ah[1] ))  - 1./(fit[0]*fit[0])) * (1./((double) V));

 //  if (flow==0 && global_c==1){ printf("%f   %f  %.15g  \n",fit[0],ct[0],(exp(-fit[0]*flow)+exp(-fit[0]*(L0-flow)))); global_c=2;}
     free(masses[0]);free(masses[1]);free(fit);free(masses);
     
}
else if(call_f==6) {
    double **masses,*fit;
    
    masses=(double**) malloc(sizeof(double*)*2);
    masses[0]=(double*) malloc(sizeof(double)*L0);
    masses[1]=(double*) malloc(sizeof(double)*L0);
    
    
    for(i=0;i<L0/2;i++){
        if(i!=0){
            ct[0]=(ah[0+i*2]+ah[0+(L0-i)*2])/2.;
            ct[1]=(ah[1+i*2]+ah[1+(L0-i)*2])/2.;
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=(ah[0+(i+1)*2]+ah[0+(L0-(i+1))*2])/2.;
            ctp[1]=(ah[1+(i+1)*2]+ah[1+(L0-(i+1))*2])/2.;
            ctm[0]=(ah[0+(i-1)*2]+ah[0+(L0-(i-1))*2])/2.;
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];

        }
        else {
            ct[0]=ah[0+i*2];
            ct[1]=ah[1+i*2];
            ct[0]=ct[0]-ah[1]*ah[1];

            ctp[0]=ah[0+(i+1)*2];
            ctp[1]=ah[1+(i+1)*2];
            ctm[0]=ah[0+(i-1)*2];
            ctp[0]=ctp[0]-ah[1]*ah[1];
            ctm[0]=ctm[0]-ah[1]*ah[1];
        }
        
        mass=log(ct[0]/ctp[0]);

        res=1;
        while(res>1e-19){
            u=1.+exp(-mass*(L0-2*i-2));
            d=1.+exp(-mass*(L0-2*i));
            tmp_mass=log( (ct[0]/ctp[0]) * (u/d)) ;
            res=tmp_mass - mass;
            mass=tmp_mass;
        }
        
        
        masses[0][i]=mass;
        masses[1][i]=effective_mass[global_k][1][i];
    }
    
     fit=constant_fit(fit_min,fit_max,masses);
  
     
     
     r[0]=ah[1]*ah[1]/(fit[0]*fit[0]*( ah[1+1*2]-((double)V)* ah[1]*ah[1] ) );
     r[0]=fit[0]*fit[0]/(2*r[0]);     

     free(masses[0]);free(masses[1]);free(fit);free(masses);
     
}
return r;
}


double **barf( int var, int order, int rep,int nconf,int flow, double *bba,double **ga)
{
double **r,*tmp,*tmp1;
double *h,*ah;
int i,j,N,alpha;

N=rep*nconf;
alpha=(order+1)*var;
ah=(double*) malloc(alpha*sizeof(double));
h=(double*) malloc(alpha*sizeof(double));
r=(double**) malloc(alpha*sizeof(double*));
for(i=0;i<alpha;i++)
{
  r[i]=(double*) calloc(order+1,sizeof(double));
  h[i]=sqrt(   ga[i][i]/((double)N*4.)  );
  ah[i]=bba[i];
}

for(i=0;i<alpha;i++)
{
  ah[i]+=h[i];
  tmp1=function(var,order,flow,ah);
  ah[i]-=h[i];ah[i]-=h[i]; 
  tmp=function(var,order,flow,ah);
  //sub_pseries(order,tmp1,tmp,r[i]);
  for(j=0;j<=order;j++){
      r[i][j]=tmp1[j]-tmp[j];
      r[i][j]/=2.*h[i];
  }
  free(tmp);free(tmp1);ah[i]+=h[i];
  //scale_pseries(order,0,1./(2.*h[i]),r[i]  );
}

free(h);free(ah);
return r;
}


void mean_value(int var, int order,int rep, int nconf,int flow,double *a)
{
    int i0,i1,i2,j,alpha,imax;
    double **ab,N;
    double *Fbb,*Fb,*tmp;

    alpha=(order+1)*var;
    N=nconf*rep;
    imax=alpha*rep;    

    for(i0=0;i0<alpha;i0++)
       abb[i0]=0;
        int i,Nfit=1;
// printf("a=%f\n",a[1+38*2]); 
    ab=(double**) malloc(rep*sizeof(double*));
    for(i1=0;i1<rep;i1++)
        ab[i1]=(double*) calloc(alpha,sizeof(double));
    
    Fb=(double*) calloc(order+1,sizeof(double));
    
    for(i0=0;i0<alpha;i0++)
        for(i1=0;i1<rep;i1++)
            for(i2=0;i2<nconf;i2++)
	     ab[i1][i0]+=a[i0+i1*alpha+i2*imax];   
    //printf("nconf=%d\t replicas=%d\t alpha=%d\n",nconf,rep,alpha);
   
    for(i0=0;i0<alpha;i0++)
        for(i1=0;i1<rep;i1++)
            abb[i0]+=ab[i1][i0];

    for(i0=0;i0<alpha;i0++)
        for(i1=0;i1<rep;i1++)
            ab[i1][i0]/=(double)nconf;

    for(i0=0;i0<alpha;i0++)        
	abb[i0]/=(double) (((double) nconf)*rep); 
    
   
    Fbb=function(var,order,flow,abb);    
   /* printf("Fbb=%f\n",Fbb[0]);*/
    if(rep==1)   for(i0=0;i0<=order;i0++) obs[i0]=Fbb[i0];
    else
    {    
        for(i1=0;i1<rep;i1++)
        {
            printf(" ok untill here  replica %d\n",i1);
            tmp=function(var,order,flow,ab[i1]);
            for(j=0;j<=order;j++){
                tmp[j]/=(double)nconf;
                Fb[j]+=tmp[j];
            }
            //scale_pseries(order,0,nconf,tmp);
            //add_pseries(order,tmp,Fb,Fb);
        }
        for(j=0;j<=order;j++)
            Fb[j]/=(double) N;
        //scale_pseries(order,0,1./((double)N),Fb);
        
        for(i0=0;i0<=order;i0++)
            obs[i0]=(((double)rep)*Fbb[i0]-Fb[i0])/(((double)rep)-1.);
    }
        
    //free_dpseries(rep-1,ab);
    for(j=0;j<=order;j++){
    free(ab[j]);
    }
}

double *Gammaf( int var, int order,double **ga,double **fa)
{
int i,j,k,alpha;
double *r;
r=(double*) calloc(order+1,sizeof(double));
alpha=var*(order+1);

for(i=0;i<=order;i++)
  for(j=0;j<alpha;j++)
    for(k=0;k<alpha;k++)
      r[i]+=fa[j][i]*fa[k][i]*ga[j][k];
   
return r;
}

void windowing(int var,int order, int rep, int nconf, int flow,double *a, double *bba)
{
    double **fbba,**tmp,*g,*tau,Caa=0;
    int count=0, i,j,i1,N,alpha;
    double S=1.5;
    
    FILE  *file_tau=NULL;
    file_tau=fopen("tau_int","w+");
    if(file_tau==NULL){
	printf("unable to open analysis file"); exit(0);
    }
    
    
    
    alpha=(order+1)*var;
    
    g=(double*) calloc(order+1,sizeof(double));
    tau=(double*) calloc(order+1,sizeof(double));
    
    N=rep*nconf;
      
     for(i=0;i<=order;i++)
    {
        CbbF[i]=0;
        w[i]=-1;
    }
  //printf("%d %d %d %d %f %f \n",var,order,rep,nconf,a[1+38*2],bba[1+38*2] ) ; 
    tmp=Gamma(0,  var,  order,  rep, nconf, a,  bba);
    fbba=barf(  var,  order,  rep, nconf, flow, bba, tmp);
    gammaFbb[0]=Gammaf(var,order,tmp,fbba);
    //add_pseries(order,CbbF,gammaFbb[0],CbbF);
    for(i1=0;i1<=order;i1++)
        CbbF[i1]+=gammaFbb[0][i1];
    Caa+=tmp[0][0];
   
//printf("\nfcall=%d\n",call_f); 
//   for(i1=0;i1<80;i1++)
//  printf("%f\t",fbba[i1][0]);
//printf("\n"); 
// printf("\n%f    %f   %f\n",gammaFbb[0][0],tmp[0][0],fbba[0][0]);
    //free_dpseries(alpha-1,tmp);    
    for(i1=0;i1<alpha;i1++)
        free(tmp[i1]);

    fprintf(file_tau,"%d   \t %d %g\n",flow,0,0.5);
    
    for(i=1;i<nconf;i++)
    {   
      
        tmp=Gamma(i,  var,  order,  rep, nconf, a,  bba);
        gammaFbb[i]=Gammaf(var,order,tmp,fbba);
        if(i==1) for(j=0;j<=order;j++)
        {
            if(gammaFbb[1][j]<0){
                /*printf("there is no autocorrelation\n");*/
                w[j]=0;count++;
            }
        }  
        if(w[0]==-1)  Caa+=2*tmp[0][0];
        //free_dpseries(alpha-1,tmp);
        for(i1=0;i1<alpha;i1++)
            free(tmp[i1]);
        for(j=0;j<=order;j++)
        {
            
            if(w[j]==-1)
            {
                CbbF[j]+=2.*gammaFbb[i][j];
                taubb_intF[j]=CbbF[j]/(2.*gammaFbb[0][j]);
 		fprintf(file_tau,"%d   \t %d %g\n",flow,i,taubb_intF[j]);
        	tau[j]=0.6;
	        if(taubb_intF[j]>0.5)
                        tau[j]=S/(  log( (2.*taubb_intF[j]+1.)/(2.*taubb_intF[j]-1.)  ));
                g[j]=exp(-((double)i)/tau[j])- (tau[j]/ (sqrt((double)(i*N) ))  );
                /*printf("g=%f",g[j]);*/
                if(g[j]<0)
                {  count++;  w[j]=i; }
                
            }
            /*if(j==0) printf("gammaFbb[%d]=%0.10f\n",i,gammaFbb[i][0]);*/
             if(count==order+1) break;
        }
        free(gammaFbb[i]);
        if(count==order+1) break;
    }
    
            

    //free_dpseries(alpha-1,fbba);
    for(i1=0;i1<alpha;i1++)
        free(fbba[i1]);

    free(g);free(tau);

        for(j=0;j<=order;j++)
        {   
            if(w[j]==-1)
            {
		printf("Windowing condition order %d failed up to W = %d\n",j,nconf-1);
                w[j]=nconf-1;
	    }
        }


    for(j=0;j<=order;j++)
    {   
        
        gammaFbb[0][j]+=CbbF[j]/((double)N);
        CbbF[j]+=CbbF[j]*(2.*w[j]+1)/((double)N);
        taubb_intF[j]=CbbF[j]/(2.*gammaFbb[0][j]);
    }
    free(gammaFbb[0]);
    fprintf(file_tau,"\n");
    fclose(file_tau);
    
}


void return_answer( int var, int order ,int rep, int nconf)
{
    int i,N;
    
    
    N=rep*nconf;
    for(i=0;i<=order;i++)
    {
        dobs[i]=CbbF[i]/((double)N);
        dobs[i]=sqrt(dobs[i]);
       
        ddobs[i]=dobs[i]*sqrt((w[i]+0.5)/N);
        dtau[i]=sqrt( (w[i]+0.5-taubb_intF[i])/ ((double)N) )*2.* taubb_intF[i] ;
    }
    
}


// imagninary element

void read_phi_field (const int size, char* filename, double *in);
void write_phi_field (const int size, char* filename, double *in);
void usage();
void finalize();
void parse_args_and_init(int argc, char** argv);

// TODO: add verbosity argument
int main (int argc, char** argv) {
  int verbose = 0;
  int order=0,alpha=2*L0,conf,var=2*L0,replicas=1,j,i,k,l,NxG=1, therm=0;
  double ***tmp,mean=0,sigma=0,*tmp1,**result, *fit;
    
    parse_args_and_init(argc,argv); 
 
    int  Nfit=1;
    corr_ave=(double***) malloc(sizeof(double**)*Nfit);
    effective_mass=(double***) malloc(sizeof(double**)*Nfit);
    matrix_element=(double***) malloc(sizeof(double**)*Nfit);
    
    for ( i = 0; i <  Nfit ; i++){
        corr_ave[i]=(double**) malloc(sizeof(double*)*2);
       effective_mass[i]=(double**) malloc(sizeof(double*)*2);
        matrix_element[i]=(double**) malloc(sizeof(double*)*2);
        for ( j = 0; j <  2 ; j++){
            corr_ave[i][j]=(double*) malloc(sizeof(double)*L0);
            effective_mass[i][j]=(double*) malloc(sizeof(double)*L0);
            matrix_element[i][j]=(double*) malloc(sizeof(double)*L0);
        }
    }
    gammaFbb=(double**) malloc(sizeof(double*)*N_tot);
    CbbF=(double*) calloc(order+1,sizeof(double));
    taubb_intF=(double*) calloc(order+1,sizeof(double));
    dtau=(double*) calloc(order+1,sizeof(double));
    dobs=(double*) calloc(order+1,sizeof(double));
    ddobs=(double*) calloc(order+1,sizeof(double));
    obs=(double*) calloc(order+1,sizeof(double));
    abb=(double*) calloc(alpha,sizeof(double));
    w=(int*) calloc(order+1,sizeof(int));
    
 
        
    
    tmp=(double***) malloc(L0*sizeof(double**));
    for ( i = 0; i <  L0 ; i++){
        tmp[i]=(double**) malloc(2*sizeof(double*));
		for ( j = 0; j <  2 ; j++){
            		tmp[i][j]=(double*) malloc(sizeof(double)*N_tot);
        }
    }
 
    tmp1=(double*) calloc((N_tot/NxG)*2*L0,sizeof(double) );

    char *cor_name=(char*) malloc(sizeof(char)*500);
    int cor_type,time;
   
    FILE  *file_in=NULL;
    file_in=fopen(input_filename,"r+");
    if(file_in==NULL){
	printf("input file not found"); exit(0);
    }
    FILE  *file_in1=NULL;
    file_in1=fopen(input_filename1,"r+");
    if(file_in1==NULL){
	printf("input file not found"); exit(0);
    }
 
    for(j=0;j<N_tot;j++)
    		for(i=0;i<L0;i++)
			  tmp[i][1][j]=rand()%10;

double a,b,c,d;
    for(j=0;j<N_tot;j++){
    		for(i=0;i<L0;i++){  

			fscanf(file_in,"%d  %lf\n ",&time	,&tmp[i][0][j]);
	                if(i==0)	fscanf(file_in1,"%lf  %lf  %lf    %lf   %lf\n ",&tmp[i+1][1][j],&tmp[i+2][1][j],&tmp[i+3][1][j],&tmp[i][1][j],&d);
	//		printf("%d  %f %f\n",i,tmp[i][0][j],tmp[i][1][j]);
/*if(i==0) printf("%f\n",tmp[i][0][j]+0.837423938138333);*/
		//printf("%d  %.15g %.15g \n",time,tmp[i][0][j],tmp[i][1][j]);
			if (time!=i){
				printf("error in reading the time coordinate\nread=%d  expected=%d \n in conf %d\n",time,i,j); exit(0);
			}
		}
        
    }
   

double **jak,mean_jak=0,**function_jak;

jak=(double**) malloc(sizeof(double*)*N_tot );
function_jak=(double**) malloc(sizeof(double*)*N_tot);

for(j=0;j<N_tot;j++)
jak[j]=(double*) calloc(2*L0,sizeof(double) );
  
    for(j=0;j<N_tot;j++){
    		for(i=0;i<L0;i++){  
    			for(l=0;l<2;l++){
                     jak[N_tot-1][l+i*2]+=tmp[i][l][j];
}
   }
}
    for(j=0;j<N_tot;j++){
    		for(i=0;i<L0;i++){  
    			for(l=0;l<2;l++){
                     jak[j][l+i*2]=(jak[N_tot-1][l+i*2]-tmp[i][l][j])/( (double) (N_tot-1) );
}
   }
}




    for(j=0;j<N_tot;j++)
    		for(i=0;i<L0;i++)
    			for(l=0;l<2;l++)
			tmp1[ l+ i*2 + (j) *2*L0 ]=tmp[i][l][j];
    
    fclose(file_in);
    FILE  *file_out=NULL;
    file_out=fopen(output_filename,"w+");
    if(file_out==NULL){
	printf("unable to open output file"); exit(0);
    }
    
    var=2*L0;

    call_f=1;
    for ( k=0; k< Nfit; k++ ){
        fprintf(file_out,"#average correlator for the mu_3 value number %d\n",k);
        fprintf(file_out,"#L0  obs dobs ddobs tau_int dtau_int  \n");
        for ( i = 0; i < L0; i++){
            mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
            windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
            return_answer( var,order, replicas, N_tot/NxG);//printf("HERE  %d %d %f\n",k,i+L0,abb[2*39]);
            fprintf(file_out,"%d \t %.15f  %.15f %.15f %.15f %.15f \n",i,obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
            corr_ave[k][0][i]=obs[0];
            corr_ave[k][1][i]=dobs[0];
        }
    }    

    
    call_f=2;
    fprintf(file_out,"\n\n");
    for ( k=0; k< Nfit; k++ ){
        fprintf(file_out,"#effective mass for the mu_3 value number %d\n",k);
        fprintf(file_out,"#L0  obs dobs ddobs tau_int dtau_int  \n");
        for ( i = 0; i < L0/2; i+=1){
            mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
            windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
            return_answer( var,order, replicas, N_tot/NxG);
            fprintf(file_out,"%d \t %.15f  %.15f %.15f %.15f %.15f \n",i,obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);

            effective_mass[k][0][i]=obs[0];
            effective_mass[k][1][i]=dobs[0];
        }
        fit=constant_fit(fit_min,fit_max,effective_mass[k]);
        fprintf(file_out,"\n\n");
        fprintf(file_out,"#fit effective mass obs dobs\n%.15f  %.15f\n",fit[0],fit[1]);
        free(fit);
    }
    
call_f=5;
for(j=0;j<N_tot;j++){
function_jak[j]=function(2*L0,0,0,jak[j]);
mean_jak+=function_jak[j][0];

}
mean_jak/=(double) N_tot;

double sigma_J=0;
for(j=0;j<N_tot;j++){

sigma_J+=(function_jak[j][0]-mean_jak   )*(function_jak[j][0]-mean_jak   );

}
sigma_J*=(N_tot-1.)/((double) N_tot);
sigma_J=sqrt(sigma_J);


printf("v_R jacknife\t %.15f      %.15f  \n",mean_jak,sigma_J);
    fprintf(file_out,"\n\n");
    call_f=4;
    fprintf(file_out,"\n\n");
    for ( k=0; k< Nfit; k++ ){
        fprintf(file_out,"#Z_phi mass for the mu_3 value number %d\n",k);
        fprintf(file_out,"#L0  obs dobs ddobs tau_int dtau_int  \n");
        for ( i = 0; i < 1; i+=1){
            mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
            windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
            return_answer( var,order, replicas, N_tot/NxG);
            fprintf(file_out," %.15f  %.15f %.15f %.15f %.15f \n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
        }
    }
    call_f=5;
    fprintf(file_out,"\n\n");
    for ( k=0; k< Nfit; k++ ){
        fprintf(file_out,"#V_R^2 mass for the mu_3 value number %d\n",k);
        fprintf(file_out,"#L0  obs dobs ddobs tau_int dtau_int  \n");
        for ( i = 0; i < 1; i+=1){
            mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
            windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
            return_answer( var,order, replicas, N_tot/NxG);
            fprintf(file_out," %.15f  %.15f %.15f %.15f %.15f \n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
        }
    }
    call_f=6;
    fprintf(file_out,"\n\n");
    for ( k=0; k< Nfit; k++ ){
        fprintf(file_out,"#lambda_R mass for the mu_3 value number %d\n",k);
        fprintf(file_out,"#L0  obs dobs ddobs tau_int dtau_int  \n");
        for ( i = 0; i < 1; i+=1){
            mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
            windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
            return_answer( var,order, replicas, N_tot/NxG);
            fprintf(file_out," %.15f  %.15f %.15f %.15f %.15f \n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
        }
    }
    
    call_f=7;
    fprintf(file_out,"\n\n");
    for ( k=0; k< Nfit; k++ ){
        fprintf(file_out,"#m_R^2 mass for the mu_3 value number %d\n",k);
        fprintf(file_out,"#L0  obs dobs ddobs tau_int dtau_int  \n");
        for ( i = 0; i < 1; i+=1){
            mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
            windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
            return_answer( var,order, replicas, N_tot/NxG);
            fprintf(file_out," %.15f  %.15f %.15f %.15f %.15f \n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
        }
    }
    
    fprintf(file_out,"\n\n");
  /*
    call_f=3;
    for ( k=0; k< Nfit; k++ ){
        global_k=k;
        fprintf(file_out,"#matrix element for the mu_3 value number %d\n",k);
        fprintf(file_out,"#L0  obs dobs ddobs tau_int dtau_int  \n");
        for(i=0;i<L0/2;i++){
            mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
            windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
            return_answer( var,order, replicas, N_tot/NxG);
            fprintf(file_out,"%d \t %.15f  %.15f %.15f %.15f %.15f \n",i,obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);

            matrix_element[k][0][i]=obs[0];
            matrix_element[k][1][i]=dobs[0];

        }
        fit=constant_fit(12,20,matrix_element[k]);
            fprintf(file_out,"\n\n");
        fprintf(file_out,"#fit_matrix element obs dobs\n%.15f  %.15f\n",fit[0],fit[1]);
        free(fit);
    }
    */    
    
    /*
    fprintf(file_out,"\n\n");
    int fit_min=9,fit_max=20;
    for ( k=0; k< Nfit; k++ ){
        fit=constant_fit(fit_min,fit_max,effective_mass[k]);
        fprintf(file_out,"#fit of the effective mass the mu_3 value number %d \n",k);
        fprintf(file_out,"#Lmin Lmax  obs dobs  \n");
        printf("%d   %d  %.15f  %.15f\n",fit_min,fit_max,fit[0],fit[1]);
    }*/
    
    /*
    for ( i = 2; i < L0/2; i+=2){
        mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
        windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
        return_answer( var,order, replicas, N_tot/NxG);
        fprintf(file_out,"%d \t %.15f  %.15f %.15f %.15f %.15f \n",i,obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);

    }
    fprintf(file_out,"\n\n");
    for ( i = 3; i < L0/2; i+=2){
        mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
        windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
        return_answer( var,order, replicas, N_tot/NxG);
        fprintf(file_out,"%d \t %.15f  %.15f %.15f %.15f %.15f \n",i,obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);

    }
    fprintf(file_out,"\n\n");
    
        mean_value(var,order,replicas,N_tot/NxG,100,tmp1);
        windowing(var,order,replicas,N_tot/NxG,100,tmp1,abb);
        return_answer( var,order, replicas, N_tot/NxG);
        fprintf(file_out," \t %.15f  %.15f %.15f %.15f %.15f \n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
  */  
/*    result=(double**)   malloc(2*sizeof(double*));
    result[0]=(double*)   malloc(L0*sizeof(double));
    result[1]=(double*)   malloc(L0*sizeof(double));
    for ( k=0; k< Nfit; k++ ){
        for ( i = 0; i < L0/2; i+=1){
            mean_value(var,order,replicas,N_tot/NxG,i,tmp1);
            windowing(var,order,replicas,N_tot/NxG,i,tmp1,abb);
            return_answer( var,order, replicas, N_tot/NxG);
printf("%d \t %.15f  %.15f %.15f %.15f %.15f \n",i,obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);

           result[0][i]=obs[0];
           result[1][i]=dobs[0];
        }
    }
    
    
        fit=constant_fit(9,20,result);
        printf("fit=%.15f  %.15f\n",fit[0],fit[1]);
  */      
    fclose(file_out);
    //finalize();
    return 0;
}

void usage(){
  printf("usage FOLLOW THE ORDER: effective_mas_gauge_log -i <infile> -o <outfile> -L <int> -T <int> -N <total configuration> -i <infile2> \n");
}

void parse_args_and_init(int argc, char** argv){
  input_filename = (char*)calloc(500,sizeof(char));
  input_filename1 = (char*)calloc(500,sizeof(char));
  output_filename = (char*)calloc(500,sizeof(char));
  int c;
  
  if(argc<7){
    usage();
    exit(1);
  }         
  
  char *tmp=(char*) calloc(2,sizeof(char));
  char *tmp1=(char*) calloc(2,sizeof(char));
  tmp=argv[1];
  tmp1=(char*) "-i";
  if(strcmp(tmp,tmp1)!=0) {usage();exit(1);}
  
  input_filename = argv[2];

  tmp=argv[3];
  tmp1=(char*) "-o";
  if(strcmp(tmp,tmp1)!=0) {usage();exit(1);}
  
  output_filename = argv[4];
  
  tmp=argv[5];
  tmp1=(char*) "-L";
  if(strcmp(tmp,tmp1)!=0) {usage();exit(1);}
  
  L1=atoi(argv[6]);
  L2=L1;
  L3=L1;
  
  tmp=argv[7];
  tmp1=(char*) "-T";
  if(strcmp(tmp,tmp1)!=0) {usage();exit(1);}
  
  L0=atoi(argv[8]);
  
  tmp=argv[9];
  tmp1=(char*) "-N";
  if(strcmp(tmp,tmp1)!=0) {usage();exit(1);}
  
  N_tot=atoi(argv[10]);
  
  tmp=argv[11];
  tmp1=(char*) "-i";
  if(strcmp(tmp,tmp1)!=0) {usage();exit(1);}
  
  input_filename1 = argv[12];
 // input_filename = argv[2];
  
 // tmp=argv[3];
  //if(tmp!="-o"){usage();exit(1);}
  
  
  /*
  while ((c = getopt(argc, argv, "h?i:o:L:T:N:")) != -1) {
    switch (c) {
      case 'L':
        L1 = L2 = L3 = atoi(optarg);
        break;
      case 'T':
        L0 = atoi(optarg);
        break;
      case 'i':
        strncpy(input_filename, optarg, 500);
        break;
      case 'o':
        strncpy(output_filename, optarg, 500);
        break;
      case 'N':
        N_tot = atoi(optarg);
        break;
      case 'h':
      case '?':
      default:
        usage();
        break;
    }
  }*/
  V = L0 * L1 * L2 * L3;
  phi_field = (double*)calloc(4*V+1,sizeof(double));
}

void finalize() {
  free(input_filename);
  free(output_filename);
  free(phi_field);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void read_phi_field (const int size, char* filename, double *in) {

  FILE *f;

  f = fopen
    (filename,
     "rb");
  if (f) {
    fread (in, sizeof(double), size, f);
    printf ("\nReading phi field configuration was successful!\n");
    fclose(f);
  }
  else{
    printf ("\n\nCANNOT OPEN PHI FIELD CONFIGURATION FILE %s !\n\n",filename);
    exit (0);
  }
}

void write_phi_field (const int size, char* filename, double *in){
  FILE *f;

  f = fopen
    (filename,
     "wb");
  if (f) {
    size_t num = fwrite (in, sizeof(double), size, f);
    printf ("Writing %lu doubles to %s was successful!\n",num,filename);
    fclose(f);
  }
  else{
    printf ("\n\nCANNOT OPEN PHI FIELD CONFIGURATION FILE %s !\n\n",filename);
    exit (0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


