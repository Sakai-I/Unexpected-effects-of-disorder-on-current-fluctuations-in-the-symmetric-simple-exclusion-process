#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <vector>
using namespace std;

random_device rd;
mt19937_64 mt64(rd());

double uniform(){
    uniform_real_distribution<double> get_rand_uni_real( 0.0, 1.0 );

    return get_rand_uni_real(mt64);
}

#define L 100 //system size
#define tau0 1. //typical time
#define T 1. //temperature
#define Emin 0. //minimum of distribution
#define Emax 3. //maximum of distribution

int main(){
    FILE *fp;
    char file[100];
    int i,j;
    double X,tau;
    double taumax;

    snprintf(file,100,"uniform/tau1.txt");

    fp = fopen(file, "w");

    taumax=0.;
    for(i=0;i<L;++i){   
        X=uniform()*(Emax-Emin)+Emin;
        tau=tau0*exp(X/T);
        fprintf(fp,"%lf\n",tau);
    }

    fprintf(fp,"\n#Emin=%.1f,Emax%.1f",Emin,Emax);

    fclose(fp);
}
