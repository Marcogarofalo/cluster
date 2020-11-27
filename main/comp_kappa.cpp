#include <iostream>
#include <complex>
#include <cstring>


int main(int argc, char **argv) {

    double msq,l0;
    double k,l;
    for ( int i = 0; i < argc; i++ ) {
        if ( ( strcmp( argv[ i ], "-msq" ) == 0 )  ) {
            msq =  atof( argv[ ++i ] );
        }
        else if ( ( strcmp( argv[ i ], "-l0" ) == 0 )  ) {
            l0 =  atof( argv[ ++i ] );
        }
        else if ( ( strcmp( argv[ i ], "-h" ) == 0 ) || ( strcmp( argv[ i ], "-help" ) == 0 ) || argc!=5 ) {
            printf( "usage:  -msq  number  -l0 number \n" );
            exit( 1 );
        }
    }
    k=-8.-msq+sqrt( (8+msq)*(8+msq)+8*l0 );
    k/=(4.*l0);

    l=l0* k*k;


    printf("%.6f  %6f \n",k,l);

    return 0;
    

}
