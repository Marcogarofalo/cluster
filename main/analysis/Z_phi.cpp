#include <iostream>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <getopt.h>

using namespace std;

// global lattice data
int L0, L1, L2, L3, V;
int verbose;
char* input_filename;
char* output_filename;

double* phi_field;

double m_sigma=0.284067335;
double **gammaFbb, *CbbF,*taubb_intF,*dtau,*dobs,*ddobs,*obs,*abb;
int *w,N_tot,every;


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
double *r,p,q,m;

p=2*sin(3.1415926535/L0 );
p*=p;
q=2*sin(2*3.1415926535/L0 );
q*=q;
//printf("%f  %f\n", p,q);
r=(double*) malloc(sizeof(double)*(order+1)); 
r[0]= (ah[2]*q-ah[1]*p) / ( ah[1] -ah[2] )  ;
//printf(" %f   %f   %f \n ",ah[0],ah[1], ah[0]-ah[1]*ah[1]);
if( flow>0 ){
m=r[0];
r[0]= (ah[0]/ah[1]) * (1./((double) V) ) * ( 1/(p+r[0])) - ( 1./(r[0]*((double) V)) );

if(flow==2){
r[0]=m/(2.*r[0]);
}
if(flow==4){
r[0]= (m_sigma*m_sigma*(ah[0]-ah[3]*ah[3]*((double)V)));

}
if(flow==3){r[0]=1./(((1./ah[1] )-(1./ah[2])) * (1./(p-q ) )) ;

}

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
    for(j=0;j<=order;j++)
        free(ab[j]);
}

double *gammaf( int var, int order,double **ga,double **fa)
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
    
    alpha=(order+1)*var;
    
    g=(double*) calloc(order+1,sizeof(double));
    tau=(double*) calloc(order+1,sizeof(double));
    
    N=rep*nconf;
      
     for(i=0;i<=order;i++)
    {
        CbbF[i]=0;
        w[i]=0;
    }
    
    tmp=Gamma(0,  var,  order,  rep, nconf, a,  bba);
    fbba=barf(  var,  order,  rep, nconf, flow, bba, tmp);
    gammaFbb[0]=gammaf(var,order,tmp,fbba);
    //add_pseries(order,CbbF,gammaFbb[0],CbbF);
    for(i1=0;i1<=order;i1++)
        CbbF[i1]+=gammaFbb[0][i1];
    Caa+=tmp[0][0];
    //free_dpseries(alpha-1,tmp);    
    for(i1=0;i1<alpha;i1++)
        free(tmp[i1]);


    
    for(i=1;i<nconf;i++)
    {   
      
        tmp=Gamma(i,  var,  order,  rep, nconf, a,  bba);
        gammaFbb[i]=gammaf(var,order,tmp,fbba);
        if(i==1) for(j=0;j<=order;j++)
            if( gammaFbb[i][j]<0) {w[j]=1;
                printf("the autocorrelatioin at time 1 is negative.\n");
                
            }
        
        if(w[0]==0)  Caa+=2*tmp[0][0];
        //free_dpseries(alpha-1,tmp);
        for(i1=0;i1<alpha;i1++)
            free(tmp[i1]);
        for(j=0;j<=order;j++)
        {   
            if(w[j]==0)
            {
                CbbF[j]+=2.*gammaFbb[i][j];
                taubb_intF[j]=CbbF[j]/(2.*gammaFbb[0][j]);
        	tau[j]=0.5;
	        if(taubb_intF[j]>0.5)
	                tau[j]=S/(  log( (2.*taubb_intF[j]+1.)/(2.*taubb_intF[j]-1.)  ));
                g[j]=exp(-((double)i)/tau[j])- (tau[j]/ (sqrt((double)(i*N) ))  );
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
            if(w[j]==0)
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
static const std::complex<double> I (0.0, 1.0);

void read_phi_field (const int size, char* filename, double *in);
void write_phi_field (const int size, char* filename, double *in);
void rotate_phi_field (double *phi_in);
double get_angle (double x, double y);
void rotate_phi_field_component (double *phi_in, int ind1, int ind2, double w);
void get_phi_field_direction (double *phi_in, double *dir, int write);
void usage();
void finalize();
void parse_args_and_init(int argc, char** argv);

// TODO: add verbosity argument
int main (int argc, char** argv) {
  int verbose = 0;
  int order=0,alpha=1,conf,var=4,replicas=1,j,i, therm=0;
  double *tmp,mean=0,sigma=0;
    int counter=0;
 double m=0.0;
    parse_args_and_init(argc,argv); 

    gammaFbb=(double**) malloc(sizeof(double*)*N_tot);
    CbbF=(double*) calloc(order+1,sizeof(double));
    taubb_intF=(double*) calloc(order+1,sizeof(double));
    dtau=(double*) calloc(order+1,sizeof(double));
    dobs=(double*) calloc(order+1,sizeof(double));
    ddobs=(double*) calloc(order+1,sizeof(double));
    obs=(double*) calloc(order+1,sizeof(double));
    abb=(double*) calloc(1,sizeof(double));
    w=(int*) calloc(order+1,sizeof(int));

        
    
    tmp=(double*) calloc((N_tot/every)*var,sizeof(double));

  // reading phi field configuration
    for ( i = 0; i < (4 * V + 1); ++i)
        phi_field[i] = 0.0;

    for ( i = 0; i < V; ++i)
        phi_field[4 * i] = 1.0;

    char* file_i=(char*) malloc(500*sizeof(char));
    
double weight, m_rot;
   // for(conf=every+therm;conf<=N_tot;conf+=every){  
       // sprintf(file_i,"%s%d",input_filename,conf);
    //    printf("%s\n",file_i);
     //   read_phi_field(4 * V + 1, file_i, phi_field);
         weight = phi_field[4 * V];
FILE *in_mag;
in_mag=fopen(input_filename,"r+");
double G2z,G2p,G2q,vev;
//printf("%s\n",input_filename);
for(conf=0;conf<therm;conf+=every)
  fscanf(in_mag,"%lf\n",&m_rot);
for(conf=every+therm;conf<=N_tot;conf+=every){  

	fscanf(in_mag,"%lf  %lf  %lf  %lf %lf \n",&G2z,&G2p,&G2q,&vev,&m_rot);
//printf("%f\n",m_rot);
        counter=(conf-therm)/every-1;
   
 
	tmp[counter*var]=G2z;
        
        tmp[1+counter*var]=G2p;
	tmp[2+counter*var]=G2q;
        tmp[3+counter*var]=vev;
 
//       printf("counter=%d   m^2  %f   \tm   %f\n",counter,tmp[counter*2],tmp[counter*2+1]);
/*        for ( i = 0; i < V; i++){
            tmp[counter]=phi_field[i*4+0];
            tmp[1+i*4+counter*4*V]=phi_field[i*4+1];
            tmp[2+i*4+counter*4*V]=phi_field[i*4+2];
            tmp[3+i*4+counter*4*V]=phi_field[i*4+3];
            
        }
        mean+=tmp[0][0][counter];
*/      
    }

        mean_value(var,order,replicas,(N_tot-therm)/every,0,tmp);
        windowing(var,order,replicas,(N_tot-therm)/every,0,tmp,abb);
        return_answer( var,order, replicas,  (N_tot-therm)/every);
        //printf("windowing procedure stops at w=%d\n",w[0]);
        //printf("mean=%.16f  sigma=%.16f \tabb=%.16f \n",mean,sigma,abb[0]);
        printf("  \t%+-.16f \t%.16f  \t%.16f  \t%.16f  \t%.16f\n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
        
        mean_value(var,order,replicas,(N_tot-therm)/every,1,tmp);
        windowing(var,order,replicas,(N_tot-therm)/every,1,tmp,abb);
        return_answer( var,order, replicas,  (N_tot-therm)/every);
        //printf("windowing procedure stops at w=%d\n",w[0]);
        //printf("mean=%.16f  sigma=%.16f \tabb=%.16f \n",mean,sigma,abb[0]);
        printf("  \t%+-.16f \t%.16f  \t%.16f  \t%.16f  \t%.16f\n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);


        mean_value(var,order,replicas,(N_tot-therm)/every,2,tmp);
        windowing(var,order,replicas,(N_tot-therm)/every,2,tmp,abb);
        return_answer( var,order, replicas,  (N_tot-therm)/every);
        //printf("windowing procedure stops at w=%d\n",w[0]);
        //printf("mean=%.16f  sigma=%.16f \tabb=%.16f \n",mean,sigma,abb[0]);
        printf("  \t%+-.16f \t%.16f  \t%.16f  \t%.16f  \t%.16f\n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);

        mean_value(var,order,replicas,(N_tot-therm)/every,3,tmp);
        windowing(var,order,replicas,(N_tot-therm)/every,3,tmp,abb);
        return_answer( var,order, replicas,  (N_tot-therm)/every);
        printf("  \t%+-.16f \t%.16f  \t%.16f  \t%.16f  \t%.16f\n",obs[0],dobs[0],ddobs[0],taubb_intF[0],dtau[0]);
        

    finalize();
    return 0;
}

void usage(){
  printf("usage: rotatephi -i<infile> -o<outfile> -L<int> -T<int> -N<total configuration> -e<saved every>\n");
}

void parse_args_and_init(int argc, char** argv){
  input_filename = (char*)calloc(500,sizeof(char));
  output_filename = (char*)calloc(500,sizeof(char));
  int c;
  
  if(argc<7){
    usage();
    exit(1);
  }
  
  while ((c = getopt(argc, argv, "h?i:o:L:T:N:e:")) != -1) {
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
      case 'e':
        every = atoi(optarg);
        break;
      case 'h':
      case '?':
      default:
        usage();
        break;
    }
  }
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
   // printf ("\nReading phi field configuration was successful!\n");
    fclose(f);
  }
  else{
   // printf ("\n\nCANNOT OPEN PHI FIELD CONFIGURATION FILE %s !\n\n",filename);
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
   // printf ("Writing %u doubles to %s was successful!\n",num,filename);
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
void rotate_phi_field (double *phi_in) {

  double angle;
  double *dir = new double[4];

  get_phi_field_direction (phi_in, dir, 1);
  angle = get_angle (dir[1], dir[0]);
  rotate_phi_field_component (phi_in, 1, 0, -angle);

  get_phi_field_direction (phi_in, dir, 0);
  angle = get_angle (dir[2], dir[0]);
  rotate_phi_field_component (phi_in, 2, 0, -angle);

  get_phi_field_direction (phi_in, dir, 0);
  angle = get_angle (dir[3], dir[0]);
  rotate_phi_field_component (phi_in, 3, 0, -angle);

  get_phi_field_direction (phi_in, dir, 1);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
double get_angle (const double x, const double y) {
  double d = sqrt (x * x + y * y);
  if (d < 1E-10)
    return 0;
  double w = asin (x / d);
  if (y < 0)
    w = M_PI - w;
  return w;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void rotate_phi_field_component (double *phi_in, const int ind1, const int ind2,
			                           const double w){
  double c = cos (w);
  double s = sin (w);

  int count = 0;
  for (int i = 0; i < (int) V; i++){
    double x = phi_in[count + ind1];
    double y = phi_in[count + ind2];
    phi_in[count + ind1] = c * x + s * y;
    phi_in[count + ind2] = -s * x + c * y;
    count += 4;
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void get_phi_field_direction (double *phi_in, double *dir, int write) {

  for (int i = 0; i < 4; ++i)
    dir[i] = 0.0;
  int counter = 0;
  for (int i = 0; i < (int) V; ++i) {
    dir[0] += phi_in[counter + 0];
    dir[1] += phi_in[counter + 1];
    dir[2] += phi_in[counter + 2];
    dir[3] += phi_in[counter + 3];
    counter += 4;
  }
  for (int i = 0; i < 4; ++i)
    dir[i] /= (double) V;
 // if (write)
  //  printf ("Higgs-Field-Direction: (%e,%e,%e,%e)\n", dir[0], dir[1], dir[2], 
    //                                                                  dir[3]);
}

