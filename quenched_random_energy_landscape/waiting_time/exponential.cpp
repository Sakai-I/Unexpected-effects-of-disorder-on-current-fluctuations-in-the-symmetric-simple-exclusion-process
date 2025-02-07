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

#define T 1.5 //temperature
#define Tg 1. //glass temperature
#define L 100 //system size
#define tau0 1. //typical time

int main(){
    FILE *fp;
    char file[100];
    int i;
    double X,tau;

    snprintf(file,100,"../data/%.1f/tau1.txt",T/Tg);
    fp = fopen(file, "w");
    for(i=0;i<L;++i){   
        X=uniform();
        tau=tau0*exp(-Tg*log(X)/T);
        fprintf(fp,"%lf\n",tau);
    }
    fclose(fp);
}